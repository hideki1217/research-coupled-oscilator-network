#pragma once

#include "system.hpp"
#include <fstream>

namespace research {
template <typename State>
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
template <typename State>
auto csv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t<State>(file_name, ",");
};
template <typename State>
auto ssv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t<State>(file_name, " ");
};
template <typename State>
auto tsv_observer_t = [](const std::string& file_name) {
  return _sv_observer_t<State>(file_name, "\t");
};
}  // namespace research