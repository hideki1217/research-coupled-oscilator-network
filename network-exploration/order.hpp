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
    // |w[0]| = |w[-1]| = 1
    for (int i = 0; i < N; i++)
      system.w[i] = std::tan(PI * (double(i + 1) / (N + 1) - 0.5));
    for (int i = 0; i < N; i++) system.w[i] /= system.w[N - 1];
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
};

struct PhaseOrder {
 private:
  coarse_grained_system_t system;
  int eval_count{0};

 public:
  PhaseOrder(double T, double tol) : system(T, tol) {}
  static PhaseOrder _default() { return PhaseOrder(1000., 1e-7); }
  auto operator()(const network_t& K, std::mt19937& rng) {
    if ((eval_count++) == 0) system.set_random_state(rng);
    system.set_network(K);
    system.burn_in(0.1);
    return system.phase_order();
  }
};