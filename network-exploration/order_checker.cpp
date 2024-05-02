#include "order.hpp"

int main() {
  std::mt19937 rng(20);
  std::uniform_real_distribution unif(0., 2 * PI);

  network_t K;
  for (int i = 0; i < K.size(); i++) std::cin >> K[i];

  std::vector<double> Ts = {500., 1000., 2000., 4000, 8000};
  for (auto T : Ts) {
    auto H = Hamiltonian(0.9, T, 1e-7);
    for (auto& x : H.state) {
      x = unif(rng);
    }
    std::cout << T << ": " << H.phase_order(K) << std::endl;
  }
}