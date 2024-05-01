#pragma once

#include <iostream>

#include "const.hpp"

namespace research {

void print_matrix(const double* mtx, const char* sep = " ") {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << mtx[i * N + j] << ((j == (N - 1)) ? "" : sep);
    }
    std::cout << std::endl;
  }
}
void print_vector(const double* mtx, const char* sep = " ") {
  for (int j = 0; j < N; j++) {
    std::cout << mtx[j] << ((j == (N - 1)) ? "" : sep);
  }
  std::cout << std::endl;
}

}  // namespace research