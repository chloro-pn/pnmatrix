#include "../include/matrix_storage_cep.h"
#include "../include/matrix.h"
#include "../include/gmres_solver.h"
#include "../include/value_compare.h"

#include "../third_party/catch.hpp"

/*
 * 1    1    1     3      6
 * 0    4   -1  *  2  ==  5
 * 2   -2    1     1      1
 */

using namespace pnmatrix;

template<typename MatrixType>
static MatrixType gmres_test() {
  MatrixType m(3, 3);
  m.set_value(1, 1, 1);
  m.set_value(1, 2, 1);
  m.set_value(1, 3, 1);
  m.set_value(2, 1, 0);
  m.set_value(2, 2, 4);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 2);
  m.set_value(3, 2, -2);
  m.set_value(3, 3, 1);

  MatrixType b(3, 1);
  b.set_value(1, 1, 6);
  b.set_value(2, 1, 5);
  b.set_value(3, 1, 1);

  gmres::option op;
  op.m = 140;
  op.rm = 1e-6;
  gmres solver(op);
  auto result = solver.solve(m, b);
  return result;
}

TEST_CASE( "gmres cep test", "[calculate]") {
  auto x = gmres_test<matrix<matrix_storage_cep<double>>>();
  REQUIRE(x.get_row() == 3);
  REQUIRE(x.get_column() == 1);
  REQUIRE(value_equal(x.get_value(1, 1), 1.0));
  REQUIRE(value_equal(x.get_value(2, 1), 2.0));
  REQUIRE(value_equal(x.get_value(3, 1), 3.0));
}

#include "../include/matrix_storage_block.h"

TEST_CASE("gmres block test", "[calculate]") {
  auto x = gmres_test<matrix<matrix_storage_block<double>>>();
  REQUIRE(x.get_row() == 3);
  REQUIRE(x.get_column() == 1);
  REQUIRE(value_equal(x.get_value(1, 1), 1.0));
  REQUIRE(value_equal(x.get_value(2, 1), 2.0));
  REQUIRE(value_equal(x.get_value(3, 1), 3.0));
}
