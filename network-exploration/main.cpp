#include <BS_thread_pool.hpp>
#include <chrono>
#include <iostream>
#include <mine.hpp>

#include "order.hpp"

struct Hamiltonian {
 public:
  const double threshold;

 private:
  PhaseOrder order_evaluator;

 public:
  Hamiltonian(double threshold, PhaseOrder&& order_evaluator)
      : threshold(threshold), order_evaluator(order_evaluator) {}
  double operator()(const network_t& K, std::mt19937& rng) {
    const auto order = order_evaluator(K, rng);

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

auto parse_args(int argc, const char** argv) {
  bool performance_mode = false;
  int iteration = 1000;
  int burn_in = 100;
  double threshold = 0.99;
  int epoch_size = 100;
  std::vector<double> betas;
  std::vector<double> scales;
  int max_threads = 8;
  double initial_k = 10.;

  {
    int cur = 0;
    while (cur < argc) {
      auto name = std::string(argv[cur]);
      if (name == "--performance") {
        performance_mode = true;
        cur++;
        continue;
      }
      if (name == "--iteration") {
        iteration = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--burn-in") {
        burn_in = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--threshold") {
        threshold = std::stod(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--epoch-size") {
        epoch_size = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--max-threads") {
        max_threads = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--initial-k") {
        initial_k = std::stod(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--repricas") {
        try {
          while (++cur < argc) {
            auto s = std::string(argv[cur]);
            auto sep_pos = s.find(':');
            auto beta = std::stod(s.substr(0, sep_pos));
            auto scale = std::stod(s.substr(sep_pos + 1));
            betas.push_back(beta);
            scales.push_back(scale);
          }
        } catch (...) {
        }
        continue;
      }

      std::cout << "Invalid arguments" << std::endl;
      std::exit(1);
    }
  }

  int num_threads = std::min(max_threads, int(betas.size()));
  if (betas.size() == 0) {
    betas.push_back(1.0);
    scales.push_back(1.0);
  }

  assert(betas.size() == scales.size());
  return std::make_tuple(performance_mode, iteration, burn_in, threshold,
                         epoch_size, scales, betas, num_threads, initial_k);
}

int main(int argc, const char** argv) {
  const auto [performance_mode, iteration, burn_in, threshold, epoch_size,
              scales, betas, num_threads, initial_k] = parse_args(argc, argv);
  const double T = 1000.;
  const double tol = 1e-7;

  std::mt19937 rng(20);
  network_t initial;
  {
    for (int i = 0; i < N * N; i++) initial[i] = initial_k / N;
    for (int i = 0; i < N; i++) initial[i * N + i] = 0;
  }
  std::vector<research::Metropolice_<network_t>> mcmcs;
  for (int i = 0; i < betas.size(); i++) {
    mcmcs.emplace_back(research::Metropolice(
        Hamiltonian(threshold, PhaseOrder(T, tol)), SymFluctuate(scales[i]),
        initial, betas[i], rng()));
  }

  BS::thread_pool pool(num_threads);
  auto update = [&]() {
    static int e = 0;
    pool.submit_loop<int>(0, mcmcs.size(),
                          [&](const int l) {
                            auto& mcmc = mcmcs[l];
                            for (int j = 0; j < epoch_size; j++) {
                              mcmc.update();
                            }
                          })
        .wait();

    // swap phase
    reprica_swap(e++ % 2 == 0, mcmcs, rng);
  };

  // Burn-In
  {
    auto start = std::chrono::system_clock::now();

    for (int e = 0; e < burn_in; e++) {
      std::cout << e + 1 << " / " << burn_in << "\r" << std::flush;
      update();
    }

    auto time = std::chrono::system_clock::now() - start;
    auto time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
    if (performance_mode) {
      std::cout << "# Param" << std::endl;
      std::cout << "T: " << T << std::endl;
      std::cout << "tol: " << tol << std::endl;
      std::cout << "N: " << N << std::endl;
      std::cout << "Num of reprica: " << betas.size() << std::endl;
      std::cout << "# Time" << std::endl;
      std::cout << time_ms << " (ms) / " << burn_in * epoch_size << "(epoch) / "
                << burn_in << " (sample)" << std::endl;
      return 0;
    }
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