#include <cassert>
#include <random>
#include <iomanip>
#include <cmath>

#include "order.hpp"

void test_phase_order() {
  static_assert(N == 2);

  std::mt19937 rng(20);
  auto evaluator = PhaseOrder::_default();

  std::vector<double> k_list;
  const int M = 101;
  for (int i = 0; i < M; i++) {
    k_list.push_back((0. * (M - 1 - i) + 4. * i) / (M - 1));
  }

  auto _expected = [](double k) {
    const double omega = 1;

    const double K = k / 2;
    if (K < omega) {
      return 0.;
    } else {
      return std::cos(0.5 * std::asin(omega / K));
    }
  };

  for (auto k : k_list) {
    network_t K;
    std::fill(K.begin(), K.end(), k / N);
    for (int i = 0; i < N; i++) K[i * N + i] = 0;

    const auto actual = evaluator(K, rng);
    const auto expected = _expected(k);
    assert(std::abs(actual - expected) < 1e-2);
  }
}

int main() { test_phase_order(); }