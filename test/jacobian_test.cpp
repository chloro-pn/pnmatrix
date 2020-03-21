#include "../include/matrix_storage_cep.h"
#include "../include/matrix.h"
#include "../include/jacobian_solver.h"
#include "../include/value_compare.h"

#include "../third_party/catch.hpp"

/*
 * 8   -3    2     3      20
 * 4   11   -1  *  2  ==  33
 * 6   3    12     1      36
 */

using namespace pnmatrix;

template<typename MatrixType>
static MatrixType jacobian_test() {
  MatrixType m(3, 3);
  m.set_value(1, 1, 8);
  m.set_value(1, 2, -3);
  m.set_value(1, 3, 2);
  m.set_value(2, 1, 4);
  m.set_value(2, 2, 11);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 6);
  m.set_value(3, 2, 3);
  m.set_value(3, 3, 12);

  MatrixType b(3, 1);
  b.set_value(1, 1, 20);
  b.set_value(2, 1, 33);
  b.set_value(3, 1, 36);

  jacobian::option op;
  op.rm = 1e-6;
  jacobian solver(op);
  auto result = solver.solve(m, b);
  return result;
}

TEST_CASE( "jacobian cep test", "[calculate]") {
  auto x = jacobian_test<matrix<matrix_storage_cep<double>>>();
  REQUIRE(x.get_row() == 3);
  REQUIRE(x.get_column() == 1);
  REQUIRE(value_equal(x.get_value(1, 1), 3.0));
  REQUIRE(value_equal(x.get_value(2, 1), 2.0));
  REQUIRE(value_equal(x.get_value(3, 1), 1.0));
}

#include "../include/matrix_storage_block.h"

TEST_CASE( "jacobian block test", "[calculate]") {
  auto x = jacobian_test<matrix<matrix_storage_block<double>>>();
  REQUIRE(x.get_row() == 3);
  REQUIRE(x.get_column() == 1);
  REQUIRE(value_equal(x.get_value(1, 1), 3.0));
  REQUIRE(value_equal(x.get_value(2, 1), 2.0));
  REQUIRE(value_equal(x.get_value(3, 1), 1.0));
}
