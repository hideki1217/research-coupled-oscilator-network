#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <mine.hpp>
#include <random>

using vec_t = std::array<double, N>;
using network_t = std::array<double, N * N>;

struct coarse_grained_system_t {
 private:
  using system_t = research::system_t;
  using state_t = system_t::state_t;

 public:
  const double T;
  const double tol;

  state_t state;

 private:
  system_t system;

 public:
  coarse_grained_system_t(double T, double tol) : T(T), tol(tol) {
    for (int i = 0; i < N; i++) {
      system.w[i] = std::tan(PI * (double(i + 1) / (N + 1) - 0.5));
    }
  }

  template <typename Rng>
  void set_random_state(Rng& rng) {
    std::uniform_real_distribution unif(0., 2 * PI);
    for (int i = 0; i < N; i++) state[i] = unif(rng);
  }
  void set_network(const network_t& K) { system.K = K; }

  void burn_in(double p = 1) {
    auto stepper = boost::numeric::odeint::make_controlled<
        boost::numeric::odeint::runge_kutta_dopri5<state_t>>(tol, tol);
    boost::numeric::odeint::integrate_adaptive(stepper, system, state, 0., T * p,
                                               1.0);
  }

  auto phase_order() {
    auto observer = research::phase_order_observer_t<system_t>();
    auto stepper = boost::numeric::odeint::make_controlled<
        boost::numeric::odeint::runge_kutta_dopri5<state_t>>(tol, tol);
    boost::numeric::odeint::integrate_const(stepper, system, state, T, 2 * T,
                                            1.0, std::ref(observer));
    return observer.value();
  }

  auto freq_order() {
    state_t start;
    std::copy(state.cbegin(), state.cend(), start.begin());

    auto stepper = boost::numeric::odeint::make_controlled<
        boost::numeric::odeint::runge_kutta_dopri5<state_t>>(tol, tol);
    boost::numeric::odeint::integrate_adaptive(stepper, system, state, T, 2 * T,
                                            1.0);

    state_t &end = state;
    vec_t mean;
    for (int i=0; i<N; i++) {
      mean[i] = (end[i] - start[i]) / T;
    }
    std::sort(mean.begin(), mean.end());

    int N_omega = 1;
    int index = 0;
    while (index < mean.size()) {
      int cur = index;
      while (cur < mean.size() && std::abs(mean[index] - mean[cur]) < 1e-2) {
        cur += 1;
      }
      N_omega = std::max(N_omega, cur - index);
      index = cur;
    }
    return (N_omega - 1) / (N - 1);
  }
};

struct PhaseOrder {
 private:
  coarse_grained_system_t system;
  int count{0};
  bool init_everytime;

 public:
  PhaseOrder(double T, double tol, bool init_everytime = false) : system(T, tol), init_everytime(init_everytime) {}
  static PhaseOrder _default() { return PhaseOrder(1000., 1e-7); }
  auto operator()(const network_t& K, std::mt19937& rng) {
    if (init_everytime || (count++) == 0) system.set_random_state(rng);
    system.set_network(K);
    system.burn_in();
    return system.phase_order();
  }
};

struct FreqOrder {
 private:
  coarse_grained_system_t system;
  int count{0};
  bool init_everytime;

 public:
  FreqOrder(double T, double tol, bool init_everytime = false)
      : system(T, tol), init_everytime(init_everytime) {}
  static FreqOrder _default() { return FreqOrder(1000., 1e-7); }
  auto operator()(const network_t& K, std::mt19937& rng) {
    if (init_everytime || (count++) == 0) system.set_random_state(rng);
    system.set_network(K);
    system.burn_in();
    return system.freq_order();
  }
};