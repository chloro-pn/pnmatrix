#include "../include/matrix_storage_block.h"
#include "../include/matrix.h"
#include "../include/gauss_seidel_solver.h"
#include "../include/value_compare.h"

#include "../third_party/catch.hpp"

/*
 * 8   -3    2     3      20
 * 4   11   -1  *  2  ==  33
 * 6   3    12     1      36
 */

using namespace pnmatrix;

static matrix<matrix_storage_block<double>> jacobian_test() {
  matrix<matrix_storage_block<double>> m(3, 3);
  m.set_value(1, 1, 8);
  m.set_value(1, 2, -3);
  m.set_value(1, 3, 2);
  m.set_value(2, 1, 4);
  m.set_value(2, 2, 11);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 6);
  m.set_value(3, 2, 3);
  m.set_value(3, 3, 12);

  matrix<matrix_storage_block<double>> b(3, 1);
  b.set_value(1, 1, 20);
  b.set_value(2, 1, 33);
  b.set_value(3, 1, 36);

  gauss_seidel::option op;
  op.rm = 1e-6;
  gauss_seidel solver(op);
  auto result = solver.solve(m, b);
  return result;
}

TEST_CASE( "guass seidel test", "[calculate]") {
  auto x = jacobian_test();
  REQUIRE(x.get_row() == 3);
  REQUIRE(x.get_column() == 1);
  REQUIRE(value_equal(x.get_value(1, 1), 3.0));
  REQUIRE(value_equal(x.get_value(2, 1), 2.0));
  REQUIRE(value_equal(x.get_value(3, 1), 1.0));
}
