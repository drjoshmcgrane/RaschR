#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

void require_finite_vector(const NumericVector& x, const char* what,
                           bool positive = false) {
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i]) || (positive && x[i] <= 0.0))
      stop("%s must contain %sfinite values", what,
           positive ? "positive " : "");
  }
}

void require_finite_matrix(const NumericMatrix& x, const char* what) {
  for (R_xlen_t i = 0; i < x.size(); ++i)
    if (!R_finite(x[i])) stop("%s must contain finite values", what);
}

NumericMatrix efrm_likelihood_impl(const NumericVector& u,
                                   const LogicalMatrix& obs,
                                   const NumericVector& score,
                                   const List& taus,
                                   const NumericVector& discs) {
  const int n = obs.nrow();
  const int J = obs.ncol();
  const int K = u.size();
  if (n < 1 || J < 1 || K < 1 || score.size() != n ||
      taus.size() != J || discs.size() != J)
    stop("incompatible EFRM likelihood inputs");
  require_finite_vector(u, "EFRM latent grid");
  require_finite_vector(score, "EFRM weighted scores");
  require_finite_vector(discs, "EFRM discriminations", true);
  for (int j = 0; j < J; ++j) {
    NumericVector tt = taus[j];
    if (tt.size() < 1) stop("every EFRM item needs at least one threshold");
    require_finite_vector(tt, "EFRM thresholds");
  }

  NumericMatrix out(n, K);
  for (int i = 0; i < n; ++i) {
    for (int q = 0; q < K; ++q) {
      out(i, q) = score[i] * u[q];
      if (!R_finite(out(i, q)))
        stop("EFRM likelihood inputs are outside the numerically representable range");
    }
  }

  std::unordered_map<std::string, std::vector<double> > cache;
  cache.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    std::string key;
    key.reserve(static_cast<std::size_t>(J));
    for (int j = 0; j < J; ++j)
      key.push_back(obs(i, j) == TRUE ? '1' : '0');

    auto found = cache.find(key);
    if (found == cache.end()) {
      std::vector<double> den(static_cast<std::size_t>(K), 0.0);
      for (int j = 0; j < J; ++j) {
        if (obs(i, j) != TRUE) continue;
        NumericVector tt = taus[j];
        const int m = tt.size();
        std::vector<double> cumulative(static_cast<std::size_t>(m + 1), 0.0);
        for (int x = 1; x <= m; ++x) {
          cumulative[static_cast<std::size_t>(x)] =
            cumulative[static_cast<std::size_t>(x - 1)] + tt[x - 1];
          if (!R_finite(cumulative[static_cast<std::size_t>(x)]))
            stop("EFRM cumulative thresholds are outside the numerically representable range");
        }

        for (int q = 0; q < K; ++q) {
          double mx = R_NegInf;
          for (int x = 0; x <= m; ++x) {
            const double lp = discs[j] *
              (u[q] * x - cumulative[static_cast<std::size_t>(x)]);
            if (!R_finite(lp))
              stop("EFRM category logits are outside the numerically representable range");
            if (lp > mx) mx = lp;
          }
          double sum_exp = 0.0;
          for (int x = 0; x <= m; ++x) {
            const double lp = discs[j] *
              (u[q] * x - cumulative[static_cast<std::size_t>(x)]);
            sum_exp += std::exp(lp - mx);
          }
          den[static_cast<std::size_t>(q)] += mx + std::log(sum_exp);
        }
      }
      found = cache.emplace(key, std::move(den)).first;
    }
    const std::vector<double>& den = found->second;
    for (int q = 0; q < K; ++q)
      out(i, q) -= den[static_cast<std::size_t>(q)];
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
NumericMatrix efrm_likelihood_cpp(NumericVector u, LogicalMatrix obs,
                                  NumericVector score, List taus,
                                  NumericVector discs) {
  return efrm_likelihood_impl(u, obs, score, taus, discs);
}

// [[Rcpp::export]]
List efrm_fit_weights_cpp(NumericMatrix L, NumericMatrix logw,
                          IntegerVector mix_idx, NumericVector count,
                          int maxit = 100, double tol = 1e-7) {
  const int n = L.nrow();
  const int K = L.ncol();
  const int H = logw.nrow();
  if (n < 1 || K < 1 || H < 1 || logw.ncol() != K ||
      mix_idx.size() != n || count.size() != n)
    stop("incompatible EFRM mixture inputs");
  if (maxit < 1) stop("EFRM mixture maxit must be positive");
  if (!R_finite(tol) || tol < 0.0)
    stop("EFRM mixture tolerance must be finite and non-negative");
  require_finite_matrix(L, "EFRM log likelihood");
  require_finite_matrix(logw, "EFRM log masses");
  require_finite_vector(count, "EFRM pattern counts", true);
  std::vector<bool> group_seen(static_cast<std::size_t>(H), false);
  for (int i = 0; i < n; ++i) {
    if (mix_idx[i] == NA_INTEGER || mix_idx[i] < 1 || mix_idx[i] > H)
      stop("invalid EFRM mixture-group index");
    group_seen[static_cast<std::size_t>(mix_idx[i] - 1)] = true;
  }
  for (int h = 0; h < H; ++h)
    if (!group_seen[static_cast<std::size_t>(h)])
      stop("empty EFRM mixture group");
  NumericMatrix current = clone(logw);

  bool converged = false;
  double step = R_PosInf;
  double loglik = R_NegInf;
  double loglik_step = R_PosInf;
  int iterations = 0;
  const bool check_convergence = R_finite(tol) && tol > 0.0;
  const auto observed_loglik = [&]() {
    double value = 0.0;
    for (int i = 0; i < n; ++i) {
      const int h = mix_idx[i] - 1;
      double mx = R_NegInf;
      for (int q = 0; q < K; ++q)
        mx = std::max(mx, L(i, q) + current(h, q));
      double denom = 0.0;
      for (int q = 0; q < K; ++q)
        denom += std::exp(L(i, q) + current(h, q) - mx);
      value += count[i] * (mx + std::log(denom));
    }
    return value;
  };
  for (int it = 0; it < maxit; ++it) {
    NumericMatrix w(H, K);
    NumericVector group_count(H);
    for (int i = 0; i < n; ++i) {
      const int h = mix_idx[i] - 1;
      if (h < 0 || h >= H) stop("invalid EFRM mixture-group index");
      double mx = R_NegInf;
      for (int q = 0; q < K; ++q)
        mx = std::max(mx, L(i, q) + current(h, q));
      double denom = 0.0;
      for (int q = 0; q < K; ++q)
        denom += std::exp(L(i, q) + current(h, q) - mx);
      for (int q = 0; q < K; ++q) {
        const double post = std::exp(L(i, q) + current(h, q) - mx) / denom;
        w(h, q) += count[i] * post;
      }
      group_count[h] += count[i];
    }

    step = 0.0;
    for (int h = 0; h < H; ++h) {
      if (group_count[h] <= 0.0) stop("empty EFRM mixture group");
      double total = 0.0;
      for (int q = 0; q < K; ++q) {
        w(h, q) = std::max(w(h, q) / group_count[h], 1e-10);
        total += w(h, q);
      }
      for (int q = 0; q < K; ++q) {
        w(h, q) /= total;
        step = std::max(step, std::abs(w(h, q) - std::exp(current(h, q))));
      }
    }
    for (int h = 0; h < H; ++h)
      for (int q = 0; q < K; ++q)
        current(h, q) = std::log(w(h, q));

    iterations = it + 1;
    if (check_convergence) {
      const double next_loglik = observed_loglik();
      loglik_step = next_loglik - loglik;
      loglik = next_loglik;
      // Adjacent grid masses may continue to exchange weight on an almost
      // flat ridge. The observed-likelihood increment is the standard EM
      // stopping criterion and is the quantity relevant to the fitted link.
      if (R_finite(loglik_step) &&
          std::abs(loglik_step) <= tol * (1.0 + std::abs(loglik))) {
        converged = true;
        break;
      }
    }
  }
  // Coordinate-ascent calls use tol = 0 for a fixed number of EM updates.
  // Their intermediate likelihood is unused; retain the final value without
  // repeating the likelihood pass at every update.
  if (!check_convergence) loglik = observed_loglik();
  return List::create(_["logw"] = current,
                      _["converged"] = converged,
                      _["step"] = step,
                      _["loglik"] = loglik,
                      _["loglik_step"] = loglik_step,
                      _["iterations"] = iterations);
}

// [[Rcpp::export]]
double efrm_negloglik_cpp(NumericVector z, NumericVector grid,
                          LogicalMatrix obs, NumericVector score,
                          List taus, NumericVector discs,
                          NumericMatrix La, NumericMatrix logw,
                          IntegerVector mix_idx, NumericVector count) {
  if (z.size() != 2) stop("the EFRM link needs scale and origin parameters");
  require_finite_vector(z, "EFRM link parameters");
  require_finite_vector(grid, "EFRM latent grid");
  NumericVector u(grid.size());
  const double ratio = std::exp(z[0]);
  if (!R_finite(ratio) || ratio <= 0.0)
    stop("EFRM link scale is outside the numerically representable range");
  for (int q = 0; q < grid.size(); ++q) {
    u[q] = ratio * grid[q] + z[1];
    if (!R_finite(u[q]))
      stop("EFRM transformed grid is outside the numerically representable range");
  }
  NumericMatrix Lb = efrm_likelihood_impl(u, obs, score, taus, discs);
  const int n = La.nrow();
  const int K = La.ncol();
  if (n < 1 || K < 1 || Lb.nrow() != n || Lb.ncol() != K ||
      logw.ncol() != K || logw.nrow() < 1 || mix_idx.size() != n ||
      count.size() != n)
    stop("incompatible EFRM objective inputs");
  require_finite_matrix(La, "EFRM first-set likelihood");
  require_finite_matrix(logw, "EFRM log masses");
  require_finite_vector(count, "EFRM pattern counts", true);

  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    if (mix_idx[i] == NA_INTEGER || mix_idx[i] < 1 ||
        mix_idx[i] > logw.nrow())
      stop("invalid EFRM mixture-group index");
    const int h = mix_idx[i] - 1;
    double mx = R_NegInf;
    for (int q = 0; q < K; ++q)
      mx = std::max(mx, La(i, q) + Lb(i, q) + logw(h, q));
    double sum_exp = 0.0;
    for (int q = 0; q < K; ++q)
      sum_exp += std::exp(La(i, q) + Lb(i, q) + logw(h, q) - mx);
    ll += count[i] * (mx + std::log(sum_exp));
  }
  if (!R_finite(ll))
    stop("EFRM objective is outside the numerically representable range");
  return -ll;
}
