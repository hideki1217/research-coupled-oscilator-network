#pragma once

#include "system.hpp"
#include <fstream>

namespace research {
template <typename State = system_t::state_t>
struct _sv_observer_t {
public:
  using state = State;
private:
  std::ofstream fout;
  const char* sep;
public:
  _sv_observer_t(const std::string& file_name, const char* sep) : fout(file_name), sep(sep) {};
  
  void operator()(const state& x, double t) {
    fout << t;
    for (int i = 0; i < x.size(); i++) {
      fout << sep << x[i];
    }
    fout << std::endl;
  }
};
auto csv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t(file_name, ",");
};
auto ssv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t(file_name, " ");
};
auto tsv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t(file_name, "\t");
};
}  // namespace research