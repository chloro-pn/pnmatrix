#pragma once
#include "type.h"
#include "value_compare.h"
#include "qr_decomposition.h"
#include <utility>
#include <cassert>
#include <cmath>

namespace pnmatrix {
class gmres {
private:
  double rm_;
  double m_;

public:
  struct option {
    double rm = 1e-6;
    double m = 30;
  };

  gmres(option op):rm_(op.rm),m_(op.m) {

  }

  template<class MatrixType>
  MatrixType solve(MatrixType& A, MatrixType& b) {
    MatrixType x0(b.get_row(), b.get_column());
    return solveinner(A, b, x0, m_);
  }

private:
  template<class MatrixType>
  MatrixType solveinner(MatrixType& A, MatrixType& b, MatrixType& x0 ,size_type restart_m) {
    using value_type = typename MatrixType::value_type;
    MatrixType result(x0.get_row(), x0.get_column());
    while (true) {
      MatrixType r0 = b - A * x0;
      value_type beta = r0.get_vector_second_norm();
      MatrixType Vm(A.get_row(), 1);
      MatrixType v1 = r0 / r0.get_vector_second_norm();
      Vm.set_column(1, v1);

      MatrixType H(restart_m + 1, restart_m);
      for (size_type m = 1; m <= restart_m; ++m) {
        H.resize(m + 1, m);
        MatrixType vsm = Vm.get_nth_column(m);
        MatrixType wm = A * vsm;
        MatrixType wmt = tr(wm);
        for (size_type i = 1; i <= m; ++i) {
          MatrixType vsi = Vm.get_sub_matrix(1, Vm.get_row(), i, 1);
          H.set_value(i, m, wmt.get_vector_inner_product(vsi));
        }
        for (size_type i = 1; i <= m; ++i) {
          MatrixType vsi = Vm.get_sub_matrix(1, Vm.get_row(), i, 1);
          wm = wm - vsi * H.get_value(i, m);
        }
        value_type h_mplus_m = wm.get_vector_second_norm();
        H.set_value(m + 1, m, h_mplus_m);
        MatrixType e1(m, 1);
        e1.set_value(1, 1, 1);

        std::pair<MatrixType, MatrixType> qr = QR<MatrixType>(H);
        MatrixType q1 = qr.first.get_nth_column(1);
        q1.delete_row(q1.get_row());
        MatrixType R = std::move(qr.second);
        R.delete_row(R.get_row());
        MatrixType mm(R.get_row(), R.get_column());
        bool error = R.inverse_with_ert(mm);
        assert(error == true);
        value_type rm = qr.first.get_value(m + 1, 1) * beta / b.get_vector_second_norm();
        if (std::abs(rm) <= rm_) {
          result = x0 + Vm * mm * q1 * beta;
          return result;
        }
        if (m == restart_m) {
          result = x0 + Vm * mm * q1 * beta;
        }
        assert(value_equal(value_type(0), h_mplus_m) == false);
        Vm.resize(A.get_row(), m + 1);
        Vm.set_column(m + 1, wm / h_mplus_m);
      }
      x0 = std::move(result);
    }
  }
};
}

