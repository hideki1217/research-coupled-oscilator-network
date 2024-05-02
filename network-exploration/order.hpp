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

  state_t state;

 private:
  vec_t w;
  bool is_first = true;

 public:
  Hamiltonian(double threshold, double T, double tol)
      : threshold(threshold), T(T), tol(tol) {
    // |w[0]| = |w[-1]| = 1
    for (int i = 0; i < N; i++)
      w[i] = std::tan(PI * (double(i + 1) / (N + 1) - 0.5));
    for (int i = 0; i < N; i++) w[i] /= w[N - 1];
  }
  double operator()(const network_t& K, std::mt19937& rng) {
    // Only the first one is initialized by random number.
    if (is_first) {
      std::uniform_real_distribution unif(0., 2 * PI);
      for (int i = 0; i < N; i++) state[i] = unif(rng);

      is_first = false;
    }

    auto order = phase_order(K);

    if (order < threshold) {
      return std::numeric_limits<double>::infinity();
    } else {
      return std::reduce(K.cbegin(), K.cend()) / N;
    }
  }

  double phase_order(const network_t& K) {
    system_t system(&w[0], &K[0]);
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
    return observer.value();
  }
};