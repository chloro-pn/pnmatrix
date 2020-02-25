#pragma once
#include "type.h"
#include "value_compare.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <utility>

namespace pnmatrix {
template<class MatrixType>
MatrixType get_household_matrix(const MatrixType& vector) {
  assert(vector.get_column() == 1);
  bool zero_vec = true;
  using value_type = typename MatrixType::value_type;
  for (size_type i = 1; i <= vector.get_row(); ++i) {
    if (!value_equal(vector.get_value(i, 1), value_type(0))) {
      zero_vec = false;
    }
  }
  assert(zero_vec == false);

  value_type norm = vector.get_vector_second_norm();
  MatrixType tmp(vector.get_row(), vector.get_column());
  tmp.set_value(1, 1, norm);
  //tmp == au
  if (tmp == vector) {
    MatrixType w(vector.get_row(), vector.get_column());
    w.set_value(2, 1, 1);
    MatrixType H = MatrixType::get_identity_matrix(vector.get_row());
    return H - w * tr(w) * 2;
  }
  MatrixType x_au = vector - tmp;
  MatrixType w = x_au / x_au.get_vector_second_norm();
  MatrixType H = MatrixType::get_identity_matrix(vector.get_row());
  H = H - w * tr(w) * 2;
  return H;
}

template<class MatrixType>
std::pair<MatrixType, MatrixType> QR(const MatrixType& matrix) {
  if(matrix.get_row() <= matrix.get_column()) {
      fprintf(stderr,"qr decomposition needs matrix(m > n)");
      std::abort();
  }
  MatrixType R(matrix);
  MatrixType Q = MatrixType::get_identity_matrix(matrix.get_row());
  for (size_type i = 1; i <= matrix.get_column(); ++i) {
    MatrixType tmp = R.get_sub_matrix(i, matrix.get_row() - i + 1, i, 1);
    MatrixType H = get_household_matrix(tmp);
    MatrixType tmp2 = MatrixType::get_identity_matrix(matrix.get_row());
    tmp2.set_value_from_matrix(i, i, H);
    Q = tmp2 * Q;
    MatrixType A_i = R.get_sub_matrix(i, matrix.get_row() - i + 1, i, matrix.get_column() - i + 1);
    auto hai = H * A_i;
    R.set_value_from_matrix(i, i, hai);
  }
//actually is Q-t and R.
return { Q,R };
}
}
