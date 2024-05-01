#include <BS_thread_pool.hpp>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <mine.hpp>

using vec_t = std::array<double, N>;
using network_t = std::array<double, N * N>;

struct Hamiltonian {
 private:
  using system_t = research::system_t;
  using state_t = system_t::state_t;

 public:
  const double T;
  const double tol;
  const double threshold;

 private:
  vec_t w;
  state_t state;
  bool is_first = true;

 public:
  Hamiltonian(double threshold, double T, double tol)
      : threshold(threshold), T(T), tol(tol) {
    // |w[0]| = |w[-1]| = 1
    for (int i = 0; i < N; i++) w[i] = std::tan(PI * ((i + 1) / (N + 1) - 0.5));
    for (int i = 0; i < N; i++) w[i] /= w[N - 1];
  }
  double operator()(const network_t& K, std::mt19937& rng) {
    system_t system(&w[0], &K[0]);

    // Only the first one is initialized by random number.
    if (is_first) {
      std::uniform_real_distribution unif(0., 2 * PI);
      for (int i = 0; i < N; i++) state[i] = unif(rng);

      is_first = false;
    }

    // Burn-In
    {
      auto stepper = boost::numeric::odeint::make_controlled<
          boost::numeric::odeint::runge_kutta_dopri5<state_t>>(tol, tol);
      boost::numeric::odeint::integrate_adaptive(stepper, system, state, 0., T,
                                                 1.0);
    }

    // Measure
    auto observer = research::phase_order_observer_t<system_t>();
    {
      auto stepper = boost::numeric::odeint::make_controlled<
          boost::numeric::odeint::runge_kutta_dopri5<state_t>>(tol, tol);
      boost::numeric::odeint::integrate_const(stepper, system, state, T, 2 * T,
                                              1.0, std::ref(observer));
    }
    auto order = observer.value();

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

    auto rebert = [=]() mutable { state[idx] = state[xdi] = prev; };
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
  const double threshold = 0.9;
  const double T = 1000.;
  const double tol = 1e-7;
  const int num_threads = 8;
  const int period = 100;
  const int burn_in = 10;
  const int iteration = 1000;
  std::vector<double> betas;
  std::vector<double> scales;

  std::mt19937 rng(20);
  network_t initial;
  {
    for (int i = 0; i < N * N; i++) initial[i] = 5. / N;
    for (int i = 0; i < N; i++) initial[i * N + i] = 0;
  }
  std::vector<research::Metropolice_<network_t>> mcmcs;
  for (int i = 0; i < betas.size(); i++) {
    mcmcs.emplace_back(research::Metropolice(Hamiltonian(threshold, T, tol),
                                             SymFluctuate(scales[i]), initial,
                                             betas[i], rng()));
  }

  BS::thread_pool pool(num_threads);

  // Burn-In
  for (int e = 0; e < burn_in; e++) {
    pool.submit_loop<int>(0, mcmcs.size(),
                          [&](const int l) {
                            auto& mcmc = mcmcs[l];
                            for (int j = 0; j < period; j++) {
                              mcmc.update();
                            }
                          })
        .wait();

    reprica_swap(e % 2 == 0, mcmcs, rng);
  }

  // Sampling
  for (int e = 0; e < iteration; e++) {
    pool.submit_loop<int>(0, mcmcs.size(),
                          [&](const int l) {
                            auto& mcmc = mcmcs[l];
                            for (int j = 0; j < period; j++) {
                              mcmc.update();
                            }
                          })
        .wait();

    // swap phase
    reprica_swap(e % 2 == 0, mcmcs, rng);
  }
}