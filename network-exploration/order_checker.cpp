#include <iomanip>

#include "order.hpp"

int main() {
  std::mt19937 rng(20);
  std::uniform_real_distribution unif(0., 2 * PI);

  network_t K;
  std::cout << "Input network (" << N << "x" << N << ")" << std::endl << "> ";
  for (int i = 0; i < K.size(); i++) std::cin >> K[i];
  std::cout << "K=" << std::reduce(K.cbegin(), K.cend()) / N << std::endl;

  std::vector<double> Ts = {100., 1000., 10000.};
  for (auto T : Ts) {
    auto system = coarse_grained_system_t(T, 1e-7);
    system.set_random_state(rng);
    system.set_network(K);
    system.burn_in(0.1);
    
    std::vector<double> results;
    for (int i=0; i<6; i++) {
      results.push_back(system.phase_order());
    }

    double mean, var;
    {
      mean = var = 0;
      for (auto& res : results) {
        mean += res;
        var += res * res;
      }
      mean /= results.size();
      var = var / results.size() - mean * mean;
    }

    std::cout << std::scientific << std::setprecision(4);
    std::cout << T << ": ";
    for (auto res : results) {
      std::cout << res << " ";
    }
    std::cout << "> " << mean << " ± √" << var << std::endl;
  }
}