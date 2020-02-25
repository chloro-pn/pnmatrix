#pragma once

#include "type.h"
#include "value_compare.h"
#include <utility>
#include <cassert>
#include <cmath>
#include <iostream>

namespace pnmatrix {
class jacobian {
private:
  double rm_;

public:
  struct option {
    double rm = 1e-6;
  };

  jacobian(option op):rm_(op.rm) {

  }

  template<class MatrixType>
  MatrixType solve(MatrixType& coeff, MatrixType& b) {
    assert(coeff.get_column() == b.get_row() && b.get_column() == 1);
    size_type x_count = coeff.get_column();
    MatrixType x_prev(x_count, 1);
    MatrixType x_next(x_count, 1);
    size_t times = 1;
    while (true) {
      for (auto row_iter = coeff.begin(); row_iter != coeff.end(); ++row_iter) {
        double result = 0;
        double sum = 0;
        size_type row = row_iter.row_index();
        for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
          sum += (*colu_iter) * x_prev.get_value(colu_iter.column_index(), 1);
        }
        result = b.get_value(row, 1) - sum + coeff.get_value(row, row) * x_prev.get_value(row, 1);
        double t = coeff.get_value(row, row);
        result = result / t;
        if (value_equal(t, 0.0) == true) {
          x_next.set_value(row, 1, 0);
        }
        else {
          x_next.set_value(row, 1, result);
        }
      }
      MatrixType tmp = coeff * x_next;
      double max_err = max_error(tmp, b);
      if (max_err <= rm_) {
        break;
      }
      else {
        std::swap(x_prev, x_next);
        ++times;
      }
    }
    return x_next;
  }

private:
  template<class MatrixType>
  double max_error(const MatrixType& m1,const MatrixType& m2) {
    assert(m1.get_row() == m2.get_row() && m1.get_column() == m2.get_column());
    double max_err = 0;
    for (size_type row = 1; row <= m1.get_row(); ++row) {
      for (size_type colu = 1; colu <= m1.get_column(); ++colu) {
        double error = m1.get_value(row, colu) - m2.get_value(row, colu);
        error = std::abs(error);
        if (error > max_err) {
          max_err = error;
        }
      }
    }
    return max_err;
  }
};
}

