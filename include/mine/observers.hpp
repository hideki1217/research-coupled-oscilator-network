#pragma once

#include "system.hpp"
#include <fstream>

namespace research {
template <typename State = system_t::state_t>
struct ssv_observer_t {
public:
  using state = State;
private:
  std::ofstream fout;
public:
  ssv_observer_t(const std::string& FileName) : fout(FileName){};
  
  void operator()(const state& x, double t) {
    static const char* sep = " ";

    fout << t;
    for (int i = 0; i < x.size(); i++) {
      fout << sep << x[i];
    }
    fout << std::endl;
  }
};
}  // namespace research