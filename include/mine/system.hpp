#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "const.hpp"

namespace research {

struct system_t {
 public:
  using state_t = std::array<double, N>;

 public:
  std::array<double, N> w;
  std::array<double, N * N> K;

 public:
  system_t(const double* _w, const double* _K) {
    std::copy_n(_w, N, w.begin());
    std::copy_n(_K, N * N, K.begin());
  }
  system_t() {
    std::fill(w.begin(), w.end(), double(0));
    std::fill(K.begin(), K.end(), double(0));
  }

  void operator()(const state_t& x, state_t& dx, double t) {
    f(&x[0], &dx[0], t);
  }

  void d(const double* x, double* D, double t) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        D[i * N + j] = K[i * N + j] * std::cos(x[j] - x[i]);
      }
    }
    for (int i = 0; i < N; i++) D[i * N + i] = 0;
    for (int i = 0; i < N; i++) {
      double reg = 0;
      for (int j = 0; j < N; j++) reg += D[i * N + j];
      D[i * N + i] = -reg;
    }
  }

  void f(const double* x, double* dx, double t) {
    for (int i = 0; i < N; i++) {
      dx[i] = w[i];
      for (int j = 0; j < N; j++) {
        dx[i] += K[i * N + j] * std::sin(x[j] - x[i]);
      }
    }
  }
};
}  // namespace research