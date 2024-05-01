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

int main() {
  network_t initial;
  {
    for (int i = 0; i < N * N; i++) initial[i] = 5. / N;
    for (int i = 0; i < N; i++) initial[i * N + i] = 0;
  }
  auto mcmc = research::Metropolice(Hamiltonian(0.9, 1000., 1e-4),
                                    SymFluctuate(1.0), initial, 1.0);
  for (int i = 0; i < 100; i++) {
    const auto res = mcmc.update();
  }
}