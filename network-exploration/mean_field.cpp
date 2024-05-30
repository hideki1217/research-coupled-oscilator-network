#include <iomanip>

#include "order.hpp"

int main() {
  std::mt19937 rng(20);
  std::uniform_real_distribution unif(0., 2 * PI);

  const double T = 1000;
  const int M = 100;
  auto system = coarse_grained_system_t(T, 1e-7);

  std::vector<double> ks;
  for (int i = 0; i < M; i++) {
    ks.push_back((double(0) * (M - 1 - i) + double(20) * i) / (M-1));
  }
  for (auto k : ks) {
    network_t K;
    for (auto& x : K) x = k / N;
    for (int i = 0; i < N; i++) K[i * N + i] = 0;

    system.set_random_state(rng);
    system.set_network(K);
    system.burn_in();

    std::cout << std::setw(6) << std::fixed;
    std::cout << k << " ";
    for (int i = 0; i < 5; i++) {
      std::cout << system.phase_order() << " ";
    }
    std::cout << std::endl;
  }
}