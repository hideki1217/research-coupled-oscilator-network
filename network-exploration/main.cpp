#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <mine/system.hpp>
#include <mine/observers.hpp>
#include <mine/utils.hpp>

int main() {
  std::cout << "Hello Network exloration" << std::endl;

  double w[N];
  for (int i = 0; i < N; i++) w[i] = 1.;
  double K[N * N];
  for (int i = 0; i < N * N; i++) K[i] = 1. / N;
  for (int i = 0; i < N; i++) K[i * N + i] = 0;

  research::print_vector(w);
  research::print_matrix(K);

  research::system_t system(w, K);
  research::system_t::state_t state;
  for (int i=0; i<N; i++) state[i] = 0;
  {
    using state_t = research::system_t::state_t;
    const double T = 1000.;

    auto observer = research::ssv_observer_t<state_t>("test.ssv");
    auto stepper = boost::numeric::odeint::make_controlled<
        boost::numeric::odeint::runge_kutta_dopri5<state_t>>(1e-4, 1e-4);
    boost::numeric::odeint::integrate_const(stepper, system,
                                               state, 0., T, 0.5, std::ref(observer));
  }
}