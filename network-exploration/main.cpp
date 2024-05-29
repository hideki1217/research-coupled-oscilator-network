#include <BS_thread_pool.hpp>
#include <bitset>
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

struct swap_result_t {
  std::bitset<N> try_swap;
  std::bitset<N> swap_status;
};
template <typename MCMC, typename Rng>
auto reprica_swap(bool mode, std::vector<MCMC>& repricas, Rng& rng) {
  auto unif = std::uniform_real_distribution(0., 1.);

  swap_result_t res;
  for (int l = mode; l + 1 < repricas.size(); l += 2) {
    res.try_swap[l].flip();

    const auto p = std::exp((repricas[l].beta - repricas[l + 1].beta) *
                            (repricas[l].energy() - repricas[l + 1].energy()));
    if (1 < p || unif(rng) < p) {
      repricas[l].swap(repricas[l + 1]);
      res.swap_status[l].flip();
    }
  }
  return res;
}

struct param_t {
  bool performance_mode = false;
  bool init_only_first = true;
  int iteration = 1000;
  int burn_in = 100;
  double threshold = 0.99;
  int epoch_size = 100;
  std::vector<double> betas;
  std::vector<double> scales;
  int num_threads = 8;
  double initial_k = 10.;

  template <typename OS>
  void print(OS& os) const {
    os << "N: " << N << std::endl;
    os << "iteration: " << iteration << std::endl;
    os << "burn_in: " << burn_in << std::endl;
    os << "threshold: " << threshold << std::endl;
    os << "epoch_size: " << epoch_size << std::endl;
    os << "num_reprica: " << betas.size() << std::endl;

    os << "betas: " << "[";
    for (auto b : betas) os << b << " ";
    os << "]" << std::endl;

    os << "scales: " << "[";
    for (auto x : scales) os << x << " ";
    os << "]" << std::endl;

    os << "num_threads: " << num_threads << std::endl;
    os << "initial_k: " << initial_k << std::endl;

    os << "init_only_first: " << int(init_only_first) << std::endl;
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
      if (name == "--init-only-first") {
        p.init_only_first = true;
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

struct stat_t {
  std::array<int, N> _try_swap;
  std::array<int, N> _swap;
  std::vector<int> _update_acc;
  std::vector<int> _update_c;

  stat_t(int num_reprica) : _update_acc(num_reprica), _update_c(num_reprica) {
    reset();
  }

  void reset() {
    std::fill(_update_acc.begin(), _update_acc.end(), 0);
    std::fill(_update_c.begin(), _update_c.end(), 0);
    std::fill(_try_swap.begin(), _try_swap.end(), 0);
    std::fill(_swap.begin(), _swap.end(), 0);
  }

  void regist_swap(const swap_result_t& result) {
    for (int i = 0; i < N; i++) {
      if (result.try_swap[i]) {
        _try_swap[i] += 1;
        _swap[i] += result.swap_status[i];
      }
    }
  }

  void regist_update(int index, bool is_accepted) {
    _update_c[index] += 1;
    _update_acc[index] += is_accepted;
  }

  double swap_rate(int index) {
    assert(index < N);
    return double(_swap[index]) / _try_swap[index];
  }

  double accepted_rate(int index) {
    return double(_update_acc[index]) / _update_c[index];
  }
};

struct update_result_t {
  swap_result_t swap_result;
  std::vector<double> accepted_rate;
};
struct async_updater_t {
  const int nstep;
  BS::thread_pool pool;
  int _c{0};

  async_updater_t(int nstep, int num_threads)
      : nstep(nstep), pool(num_threads) {}
  template <typename Rng, typename MCMC>
  void operator()(std::vector<MCMC>& repricas, Rng& rng, stat_t& stat) {
    pool.submit_loop<int>(0, repricas.size(),
                          [&](const int l) {
                            auto& mcmc = repricas[l];
                            for (int j = 0; j < nstep; j++) {
                              const auto status = mcmc.update();
                              stat.regist_update(
                                  l, status == MCMC::Result::Accepted);
                            }
                          })
        .wait();

    // swap phase
    const auto swap_res = reprica_swap(_c++ % 2 == 0, repricas, rng);
    stat.regist_swap(swap_res);
  }

  template <typename Rng, typename MCMC>
  void operator()(std::vector<MCMC>& repricas, Rng& rng) {
    pool.submit_loop<int>(0, repricas.size(),
                          [&](const int l) {
                            auto& mcmc = repricas[l];
                            for (int j = 0; j < nstep; j++) {
                              mcmc.update();
                            }
                          })
        .wait();

    // swap phase
    reprica_swap(_c++ % 2 == 0, repricas, rng);
  }
};

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
        Hamiltonian(p.threshold, PhaseOrder(T, tol, p.init_only_first)), SymFluctuate(p.scales[i]),
        initial, p.betas[i], rng()));
  }

  async_updater_t updater(p.epoch_size, p.num_threads);
  stat_t stat(mcmcs.size());

  {
    std::ofstream param_f("param.yaml");
    p.print(param_f);
  }

  // Burn-In
  {
    auto start = std::chrono::system_clock::now();

    for (int e = 0; e < p.burn_in; e++) {
      updater(mcmcs, rng);
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
      std::cout << time_ms << " (ms) / " << p.burn_in * p.epoch_size
                << "(epoch) / " << p.burn_in << " (sample)" << std::endl;
      return 0;
    }
  }

  try {
    std::ofstream network_f("network.ssv");
    // Sampling
    for (int e = 0; e < p.iteration; e++) {
      updater(mcmcs, rng, stat);

      for (auto& mcmc : mcmcs) {
        for (int i = 0; i < mcmc.state().size(); i++) {
          network_f << mcmc.state()[i]
                    << ((i != mcmc.state().size() - 1) ? " " : "");
        }
        network_f << std::endl;
      }

      if (e % 100 == 0) {
        std::cout << "e=" << e << std::endl;

        std::cout << stat.swap_rate(0);
        for (int i = 1; i < N; i++) {
          std::cout << " " << stat.swap_rate(i);
        }
        std::cout << std::endl;

        std::cout << stat.accepted_rate(0);
        for (int i = 1; i < N; i++) {
          std::cout << " " << stat.accepted_rate(i);
        }
        std::cout << std::endl;

        stat.reset();
      }
    }
  } catch (...) {
  }
}