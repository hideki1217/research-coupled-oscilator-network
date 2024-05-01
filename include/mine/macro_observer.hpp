#include <cassert>
#include <cmath>

#include "system.hpp"

namespace research {
template<typename System>
struct phase_order_observer_t {
 public:
  using target_t = System;
  using state_t = typename System::state_t;

 private:
  double m_cos;
  double m_sin;
  int count;

 public:
  phase_order_observer_t() : m_cos(0), m_sin(0), count(0) {}

 public:
  void operator()(const state_t& x, double t) {
    double cos = 0;
    double sin = 0;
    {
      for (int i = 0; i < N; i++) {
        cos += std::cos(x[i]);
        sin += std::sin(x[i]);
      }
      cos /= N;
      sin /= N;
    }

    m_cos += (cos - m_cos) / (count + 1);
    m_sin += (sin - m_sin) / (count + 1);
    count += 1;
  }

  auto value() {
    const auto res = std::sqrt(m_cos * m_cos + m_sin * m_sin);
    if (std::isnormal(res)) {
      return res;
    } else {
      assert(std::isfinite(res));
      return 0.;
    }
  }
};

}  // namespace research