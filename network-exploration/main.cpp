#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <mine.hpp>

int main() {
  std::cout << "Hello Network exloration" << std::endl;

  double w[N];
  for (int i = 0; i < N; i++) w[i] = 1.;
  double K[N * N];
  for (int i = 0; i < N * N; i++) K[i] = 1. / N;
  for (int i = 0; i < N; i++) K[i * N + i] = 0;

  research::system_t system(w, K);
  research::system_t::state_t state;
  for (int i = 0; i < N; i++) state[i] = 0;
  
  std::vector<double> Ts;
  for (int i=0; i<30; i++) Ts.push_back(1 * std::pow(1.5, i));
  for (auto T: Ts) 
  {
    using state_t = research::system_t::state_t;

    auto observer = research::phase_order_observer_t();
    auto stepper = boost::numeric::odeint::make_controlled<
        boost::numeric::odeint::runge_kutta_dopri5<state_t>>(1e-4, 1e-4);
    boost::numeric::odeint::integrate_const(stepper, system, state, 0., T, 1.0,
                                            std::ref(observer));
    std::cout << T << " " << observer.value() << " " << T * observer.value() << std::endl;
  }
}