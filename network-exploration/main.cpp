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

struct AsymFluctuate {
 public:
  const double scale;

 private:
  std::vector<std::pair<int, int>> index_list;

 public:
  AsymFluctuate(double scale) : scale(scale) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (i != j) index_list.push_back({i, j});
      }
    }
  }
  std::function<void()> operator()(network_t& state, std::mt19937& rng) {
    std::normal_distribution normal(0., scale);
    std::uniform_int_distribution indexer(0, int(index_list.size() - 1));

    const auto [i, j] = index_list[indexer(rng)];
    const auto noise = normal(rng);

    const auto idx = i * N + j;
    const auto prev = state[idx];

    state[idx] += noise;

    auto rebert = [=, &state]() mutable { state[idx] = prev; };
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

struct param_t {
  bool performance_mode = false;
  int iteration = 1000;
  int burn_in = 100;
  double threshold = 0.99;
  int epoch_size = 100;
  std::vector<double> betas;
  std::vector<double> scales;
  int num_threads = 8;
  double initial_k = 10.;

  void print() const {
    std::cout << "# " << "iteration: " << iteration << std::endl;
    std::cout << "# " << "burn_in: " << burn_in << std::endl;
    std::cout << "# " << "threshold: " << threshold << std::endl;
    std::cout << "# " << "epoch_size: " << epoch_size << std::endl;
    std::cout << "# " << "num_reprica: " << betas.size() << std::endl;
    
    std::cout << "# " << "betas: " << "[";
    for (auto b : betas) std::cout << b << " ";
    std::cout << "]" << std::endl;
    
    std::cout << "# " << "scales: " << "[";
    for (auto x : scales) std::cout << x << " ";
    std::cout << "]" << std::endl;

    std::cout << "# " << "num_threads: " << num_threads << std::endl;
    std::cout << "# " << "initial_k: " << initial_k << std::endl;
  }
};

auto parse_args(int argc, const char** argv) {
  param_t p;

  {
    int cur = 1;
    while (cur < argc) {
      auto name = std::string(argv[cur]);
      if (name == "--performance") {
        p.performance_mode = true;
        cur++;
        continue;
      }
      if (name == "--iteration") {
        p.iteration = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--burn-in") {
        p.burn_in = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--threshold") {
        p.threshold = std::stod(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--epoch-size") {
        p.epoch_size = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--max-threads") {
        p.num_threads = std::stoi(argv[++cur]);
        cur++;
        continue;
      }
      if (name == "--initial-k") {
        p.initial_k = std::stod(argv[++cur]);
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
            p.betas.push_back(beta);
            p.scales.push_back(scale);
          }
        } catch (...) {
        }
        continue;
      }

      std::cout << "Invalid arguments: " << argv[cur] << std::endl;
      std::exit(1);
    }
  }

  p.num_threads = std::min(p.num_threads, int(p.betas.size()));
  if (p.betas.size() == 0) {
    p.betas.push_back(1.0);
    p.scales.push_back(1.0);
  }

  assert(p.betas.size() == p.scales.size());
  return p;
}

int main(int argc, const char** argv) {
  const auto p = parse_args(argc, argv);
  const double T = 1000.;
  const double tol = 1e-7;

  std::mt19937 rng(20);
  network_t initial;
  {
    for (int i = 0; i < N * N; i++) initial[i] = p.initial_k / N;
    for (int i = 0; i < N; i++) initial[i * N + i] = 0;
  }
  std::vector<research::Metropolice_<network_t>> mcmcs;
  for (int i = 0; i < p.betas.size(); i++) {
    mcmcs.emplace_back(research::Metropolice(
        Hamiltonian(p.threshold, PhaseOrder(T, tol)), SymFluctuate(p.scales[i]),
        initial, p.betas[i], rng()));
  }

  BS::thread_pool pool(p.num_threads);
  auto update = [&]() {
    static int e = 0;
    pool.submit_loop<int>(0, mcmcs.size(),
                          [&](const int l) {
                            auto& mcmc = mcmcs[l];
                            for (int j = 0; j < p.epoch_size; j++) {
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

    for (int e = 0; e < p.burn_in; e++) {
      update();
    }

    auto time = std::chrono::system_clock::now() - start;
    auto time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
    if (p.performance_mode) {
      std::cout << "# Param" << std::endl;
      std::cout << "T: " << T << std::endl;
      std::cout << "tol: " << tol << std::endl;
      std::cout << "N: " << N << std::endl;
      std::cout << "Num of reprica: " << p.betas.size() << std::endl;
      std::cout << "# Time" << std::endl;
      std::cout << time_ms << " (ms) / " << p.burn_in * p.epoch_size << "(epoch) / "
                << p.burn_in << " (sample)" << std::endl;
      return 0;
    }
  }

  p.print();

  // Sampling
  for (int e = 0; e < p.iteration; e++) {
    update();

    for (auto& mcmc : mcmcs) {
      for (auto& K_ij : mcmc.state()) {
        std::cout << K_ij << " ";
      }
      std::cout << std::endl;
    }
  }
}