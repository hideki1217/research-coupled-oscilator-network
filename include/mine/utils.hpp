#pragma once

#include <iostream>

#include "const.hpp"

namespace research {

void print_matrix(int n, const double* mtx, const char* sep = " ") {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << mtx[i * n + j] << ((j == (n - 1)) ? "" : sep);
    }
    std::cout << std::endl;
  }
}
void print_vector(int n, const double* mtx, const char* sep = " ") {
  for (int j = 0; j < n; j++) {
    std::cout << mtx[j] << ((j == (n - 1)) ? "" : sep);
  }
  std::cout << std::endl;
}

}  // namespace research