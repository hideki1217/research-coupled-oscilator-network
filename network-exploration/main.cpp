#include <BS_thread_pool.hpp>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <mine.hpp>

#include "order.hpp"

struct Hamiltonian {
 public:
  const double threshold;
  coarse_grained_system_t system;

 private:
  int eval_count{0};

 public:
  Hamiltonian(double threshold, double T, double tol)
      : threshold(threshold), system(T, tol) {}
  double operator()(const network_t& K, std::mt19937& rng) {
    if ((eval_count++) == 0) system.set_random_state(rng);
    system.set_network(K);
    system.burn_in();
    auto order = system.phase_order();

    if (order < threshold) {
      return std::numeric_limits<double>::infinity();
    } else {
      return std::reduce(K.cbegin(), K.cend()) / N;
    }
  }
};

struct SymFluctuate {
 public:
  const double scale;

 private:
  std::vector<std::pair<int, int>> index_list;

 public:
  SymFluctuate(double scale) : scale(scale) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (i > j) index_list.push_back({i, j});
      }
    }
  }
  std::function<void()> operator()(network_t& state, std::mt19937& rng) {
    std::normal_distribution normal(0., scale);
    std::uniform_int_distribution indexer(0, int(index_list.size() - 1));

    const auto [i, j] = index_list[indexer(rng)];
    const auto noise = normal(rng);

    const auto idx = i * N + j;
    const auto xdi = j * N + i;
    const auto prev = state[idx];

    state[idx] += noise;
    state[xdi] += noise;

    auto rebert = [=, &state]() mutable { state[idx] = state[xdi] = prev; };
    return rebert;
  }
};

template <typename MCMC, typename Rng>
void reprica_swap(bool mode, std::vector<MCMC>& repricas, Rng& rng) {
  auto unif = std::uniform_real_distribution(0., 1.);

  for (int l = mode; l + 1 < repricas.size(); l += 2) {
    const auto p = std::exp((repricas[l].beta - repricas[l + 1].beta) *
                            (repricas[l].energy() - repricas[l + 1].energy()));
    if (1 < p || unif(rng) < p) {
      repricas[l].swap(repricas[l + 1]);
    }
  }
}

int main() {
  const double threshold = 0.99;
  const double T = 1000.;
  const double tol = 1e-7;
  const int num_threads = 8;
  const int period = 100;
  const int burn_in = 100;
  const int iteration = 1000;
  std::vector<double> betas = {1.0, 1.5};
  std::vector<double> scales = {1.0, 1 / 1.5};

  std::mt19937 rng(20);
  network_t initial;
  {
    for (int i = 0; i < N * N; i++) initial[i] = 10. / N;
    for (int i = 0; i < N; i++) initial[i * N + i] = 0;
  }
  std::vector<research::Metropolice_<network_t>> mcmcs;
  for (int i = 0; i < betas.size(); i++) {
    mcmcs.emplace_back(research::Metropolice(Hamiltonian(threshold, T, tol),
                                             SymFluctuate(scales[i]), initial,
                                             betas[i], rng()));
  }

  BS::thread_pool pool(num_threads);
  auto update = [&]() {
    static int e = 0;
    pool.submit_loop<int>(0, mcmcs.size(),
                          [&](const int l) {
                            auto& mcmc = mcmcs[l];
                            for (int j = 0; j < period; j++) {
                              mcmc.update();
                            }
                          })
        .wait();

    // swap phase
    reprica_swap(e++ % 2 == 0, mcmcs, rng);
  };

  // Burn-In
  for (int e = 0; e < burn_in; e++) {
    update();
  }

  // Sampling
  for (int e = 0; e < iteration; e++) {
    update();

    for (auto& mcmc : mcmcs) {
      for (auto& K_ij : mcmc.state()) {
        std::cout << K_ij << " ";
      }
    }
    std::cout << std::endl;
  }
}