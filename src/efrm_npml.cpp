#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

NumericMatrix efrm_likelihood_impl(const NumericVector& u,
                                   const LogicalMatrix& obs,
                                   const NumericVector& score,
                                   const List& taus,
                                   const NumericVector& discs) {
  const int n = obs.nrow();
  const int J = obs.ncol();
  const int K = u.size();
  if (score.size() != n || taus.size() != J || discs.size() != J)
    stop("incompatible EFRM likelihood inputs");

  NumericMatrix out(n, K);
  for (int i = 0; i < n; ++i)
    for (int q = 0; q < K; ++q)
      out(i, q) = score[i] * u[q];

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
        for (int x = 1; x <= m; ++x)
          cumulative[static_cast<std::size_t>(x)] =
            cumulative[static_cast<std::size_t>(x - 1)] + tt[x - 1];

        for (int q = 0; q < K; ++q) {
          double mx = R_NegInf;
          for (int x = 0; x <= m; ++x) {
            const double lp = discs[j] *
              (u[q] * x - cumulative[static_cast<std::size_t>(x)]);
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
  if (logw.ncol() != K || mix_idx.size() != n || count.size() != n)
    stop("incompatible EFRM mixture inputs");
  NumericMatrix current = clone(logw);

  bool converged = false;
  double step = R_PosInf;
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
    if (step < tol) {
      converged = true;
      break;
    }
  }
  return List::create(_["logw"] = current,
                      _["converged"] = converged,
                      _["step"] = step);
}

// [[Rcpp::export]]
double efrm_negloglik_cpp(NumericVector z, NumericVector grid,
                          LogicalMatrix obs, NumericVector score,
                          List taus, NumericVector discs,
                          NumericMatrix La, NumericMatrix logw,
                          IntegerVector mix_idx, NumericVector count) {
  if (z.size() != 2) stop("the EFRM link needs scale and origin parameters");
  NumericVector u(grid.size());
  const double ratio = std::exp(z[0]);
  for (int q = 0; q < grid.size(); ++q)
    u[q] = ratio * grid[q] + z[1];
  NumericMatrix Lb = efrm_likelihood_impl(u, obs, score, taus, discs);
  const int n = La.nrow();
  const int K = La.ncol();
  if (Lb.nrow() != n || Lb.ncol() != K || mix_idx.size() != n ||
      count.size() != n)
    stop("incompatible EFRM objective inputs");

  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    const int h = mix_idx[i] - 1;
    if (h < 0 || h >= logw.nrow())
      stop("invalid EFRM mixture-group index");
    double mx = R_NegInf;
    for (int q = 0; q < K; ++q)
      mx = std::max(mx, La(i, q) + Lb(i, q) + logw(h, q));
    double sum_exp = 0.0;
    for (int q = 0; q < K; ++q)
      sum_exp += std::exp(La(i, q) + Lb(i, q) + logw(h, q) - mx);
    ll += count[i] * (mx + std::log(sum_exp));
  }
  return -ll;
}
