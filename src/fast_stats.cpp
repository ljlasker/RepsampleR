#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

static NumericVector insert_sorted_impl(const NumericVector& sorted_x, double value) {
  const int n = sorted_x.size();
  NumericVector out(n + 1);

  int pos = std::lower_bound(sorted_x.begin(), sorted_x.end(), value) - sorted_x.begin();

  for (int i = 0; i < pos; ++i) {
    out[i] = sorted_x[i];
  }
  out[pos] = value;
  for (int i = pos; i < n; ++i) {
    out[i + 1] = sorted_x[i];
  }

  return out;
}

// [[Rcpp::export]]
NumericVector insert_sorted_cpp(NumericVector sorted_x, double value) {
  return insert_sorted_impl(sorted_x, value);
}

// [[Rcpp::export]]
double ks2_d_with_candidate_cpp(NumericVector sample_sorted,
                                double candidate,
                                NumericVector pop_sorted) {
  NumericVector sample_plus = insert_sorted_impl(sample_sorted, candidate);

  const int n1 = sample_plus.size();
  const int n2 = pop_sorted.size();

  int i = 0;
  int j = 0;
  double d = 0.0;

  while (i < n1 || j < n2) {
    double v;

    if (j >= n2 || (i < n1 && sample_plus[i] < pop_sorted[j])) {
      v = sample_plus[i];
    } else if (i >= n1 || pop_sorted[j] < sample_plus[i]) {
      v = pop_sorted[j];
    } else {
      v = sample_plus[i];
    }

    while (i < n1 && sample_plus[i] == v) {
      ++i;
    }
    while (j < n2 && pop_sorted[j] == v) {
      ++j;
    }

    const double cdf1 = static_cast<double>(i) / static_cast<double>(n1);
    const double cdf2 = static_cast<double>(j) / static_cast<double>(n2);
    const double diff = std::fabs(cdf1 - cdf2);

    if (diff > d) {
      d = diff;
    }
  }

  return d;
}

// [[Rcpp::export]]
double ks1_norm_d_with_candidate_cpp(NumericVector sample_sorted,
                                     double candidate,
                                     double mean,
                                     double sd) {
  NumericVector sample_plus = insert_sorted_impl(sample_sorted, candidate);

  const int n = sample_plus.size();
  double d = 0.0;

  for (int i = 0; i < n; ++i) {
    const double f = R::pnorm(sample_plus[i], mean, sd, 1, 0);
    const double d_plus = static_cast<double>(i + 1) / static_cast<double>(n) - f;
    const double d_minus = f - static_cast<double>(i) / static_cast<double>(n);

    if (d_plus > d) {
      d = d_plus;
    }
    if (d_minus > d) {
      d = d_minus;
    }
  }

  return d;
}

// [[Rcpp::export]]
double ks1_cdf_d_with_candidate_cpp(NumericVector sample_sorted,
                                    NumericVector sample_cdf_sorted,
                                    double candidate,
                                    double candidate_cdf) {
  const int n0 = sample_sorted.size();
  if (sample_cdf_sorted.size() != n0) {
    stop("`sample_cdf_sorted` must have the same length as `sample_sorted`.");
  }

  NumericVector sample_plus_cdf(n0 + 1);
  const int pos = std::lower_bound(sample_sorted.begin(), sample_sorted.end(), candidate) - sample_sorted.begin();

  for (int i = 0; i < pos; ++i) {
    sample_plus_cdf[i] = sample_cdf_sorted[i];
  }
  sample_plus_cdf[pos] = candidate_cdf;
  for (int i = pos; i < n0; ++i) {
    sample_plus_cdf[i + 1] = sample_cdf_sorted[i];
  }

  const int n = sample_plus_cdf.size();
  double d = 0.0;

  for (int i = 0; i < n; ++i) {
    const double f = sample_plus_cdf[i];
    const double d_plus = static_cast<double>(i + 1) / static_cast<double>(n) - f;
    const double d_minus = f - static_cast<double>(i) / static_cast<double>(n);

    if (d_plus > d) {
      d = d_plus;
    }
    if (d_minus > d) {
      d = d_minus;
    }
  }

  return d;
}
