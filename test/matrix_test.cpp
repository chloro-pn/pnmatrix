#include "../third_party/catch.hpp"
#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/value_compare.h"
#include "../include/matrix_storage_block.h"

using namespace pnmatrix;

TEST_CASE("matrix set and get value test ","[matrix]") {
  SECTION("cep") {
    matrix<matrix_storage_cep<double>> m1(3, 3);
    REQUIRE(value_equal(m1.get_value(1, 1), 0.0));
    m1.set_value(2, 1, 1.05);
    REQUIRE(value_equal(m1.get_value(2, 1), 1.05));
    m1.add_value(2, 1, 1.1);
    REQUIRE(value_equal(m1.get_value(2, 1), 1.05 + 1.1));
  }
  SECTION("block") {
    matrix<matrix_storage_block<double>> m1(3, 3);
    REQUIRE(value_equal(m1.get_value(1, 1), 0.0));
    m1.set_value(2, 1, 1.05);
    REQUIRE(value_equal(m1.get_value(2, 1), 1.05));
    m1.add_value(2, 1, 1.1);
    REQUIRE(value_equal(m1.get_value(2, 1), 1.05 + 1.1));
  }
}

TEST_CASE("matrix iterator test", "[matrix]") {
  matrix<matrix_storage_cep<double>> m1(4, 4);
  std::vector<double> nodes;
  m1.set_value(1, 1, 1.01);
  nodes.push_back(1.01);
  m1.set_value(1, 4, 2.12);
  nodes.push_back(2.12);
  m1.set_value(2, 1, 2.2);
  nodes.push_back(2.2);
  m1.set_value(2, 3, 0.98);
  nodes.push_back(0.98);
  m1.set_value(3, 3, 3.1);
  nodes.push_back(3.1);
  m1.set_value(4, 2, 1.04);
  nodes.push_back(1.04);

  auto iter = nodes.begin();
  for (auto row = m1.begin(); row != m1.end(); ++row) {
    for (auto col = row.begin(); col != row.end(); ++col) {
      REQUIRE(value_equal(*col, *iter));
      ++iter;
    }
  }
}

TEST_CASE("matrix copy and move test", "[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 3);
  auto m2 = m1;
  bool e = (m1 == m2);
  REQUIRE(e == true);
  m1.set_value(1, 1, 2.12);
  e = (m1 == m2);
  REQUIRE(e == false);
  m2 = m1;
  e = (m1 == m2);
  REQUIRE(e == true);
  auto m3(m1);
  e = (m1 == m3);
  REQUIRE(e == true);
  m3.set_value(2, 2, 1.04);
  m2 = std::move(m3);
  REQUIRE(value_equal(m2.get_value(2, 2), 1.04));
}

TEST_CASE("matrix delete row and column test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 3);
  m1.set_value(1, 1, 2.14);
  m1.set_value(2, 1, 3.3);
  m1.set_value(3, 3, 1.08);

  SECTION("delete last row and column") {
    m1.delete_row(3);
    REQUIRE(m1.get_row() == 2);
    REQUIRE(value_equal(m1.get_value(2, 1), 3.3));
    m1.delete_column(3);
    REQUIRE(m1.get_column() == 2);
    REQUIRE(value_equal(m1.get_value(1, 1), 2.14));
  }

  SECTION("delete random row and column") {
    m1.delete_row(1);
    REQUIRE(m1.get_row() == 2);
    REQUIRE(value_equal(m1.get_value(1, 1), 3.3));
    REQUIRE(value_equal(m1.get_value(2, 3), 1.08));
    m1.delete_column(2);
    REQUIRE(m1.get_column() == 2);
    REQUIRE(value_equal(m1.get_value(2, 2), 1.08));
  }
}

TEST_CASE("matrix resize test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(4, 4);
  m1.set_value(2, 3, 1.4);
  m1.set_value(3, 4, 1.1);
  m1.set_value(4, 4, 0.97);

  SECTION("row and column small") {
    m1.resize(3, 2);
    REQUIRE(m1.get_row() == 3);
    REQUIRE(m1.get_column() == 2);
    REQUIRE(value_equal(m1.get_value(3, 2), 0.0));
  }

  SECTION("row big and column small") {
    m1.resize(5, 3);
    REQUIRE(m1.get_row() == 5);
    REQUIRE(m1.get_column() == 3);
    REQUIRE(value_equal(m1.get_value(2, 3), 1.4));
    REQUIRE(value_equal(m1.get_value(5, 1), 0.0));
  }

  SECTION("row small and column big") {
    m1.resize(3, 5);
    REQUIRE(m1.get_row() == 3);
    REQUIRE(m1.get_column() == 5);
    REQUIRE(value_equal(m1.get_value(3, 4), 1.1));
    REQUIRE(value_equal(m1.get_value(3, 5), 0.0));
  }
}

TEST_CASE("matrix operator test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 3);
  m1.set_value(1, 1, 2.01);
  m1.set_value(1, 2, 4.5);
  m1.set_value(2, 2, 4.1);
  m1.set_value(3, 3, 1.4);
  SECTION("operator * value") {
    auto m2 = m1 * 2.0;
    REQUIRE(m2.get_row() == 3);
    REQUIRE(m2.get_column() == 3);
    REQUIRE(value_equal(m2.get_value(1, 1), 4.02));
    REQUIRE(value_equal(m2.get_value(2, 2), 8.2));
  }

  SECTION("operator / value") {
    auto m2 = m1 / 1.5;
    REQUIRE(m2.get_row() == 3);
    REQUIRE(m2.get_column() == 3);
    REQUIRE(value_equal(m2.get_value(1, 1), 2.01 / 1.5));
    REQUIRE(value_equal(m2.get_value(3, 3), 1.4 / 1.5));
  }

  SECTION("operator + matrix") {
    matrix<matrix_storage_cep<double>> m2(3, 3);
    m2.set_value(1, 2, 5.5);
    m2.set_value(3, 2, 1.2);
    auto m3 = m1 + m2;
    REQUIRE(m3.get_row() == 3);
    REQUIRE(m3.get_column() == 3);
    REQUIRE(value_equal(m3.get_value(1, 2), 10.0));
    REQUIRE(value_equal(m3.get_value(3, 1), 0.0));
    REQUIRE(value_equal(m3.get_value(1, 1), 2.01));
    REQUIRE(value_equal(m3.get_value(3, 2), 1.2));
  }

  SECTION("operator - matrix") {
    matrix<matrix_storage_cep<double>> m2(3, 3);
    m2.set_value(1, 2, 5.5);
    m2.set_value(3, 2, 1.2);
    auto m3 = m1 - m2;
    REQUIRE(m3.get_row() == 3);
    REQUIRE(m3.get_column() == 3);
    REQUIRE(value_equal(m3.get_value(1, 2), -1.0));
    REQUIRE(value_equal(m3.get_value(3, 1), 0.0));
    REQUIRE(value_equal(m3.get_value(1, 1), 2.01));
    REQUIRE(value_equal(m3.get_value(3, 2), -1.2));
  }

  SECTION("operator * matrix") {
    matrix<matrix_storage_cep<double>> m2(3, 1);
    m2.set_value(2, 1, 0.4);
    auto m3 = m1 * m2;
    REQUIRE(m3.get_row() == 3);
    REQUIRE(m3.get_column() == 1);
    REQUIRE(value_equal(m3.get_value(1, 1), 4.5 * 0.4));
    REQUIRE(value_equal(m3.get_value(2, 1), 4.1 * 0.4));
  }
}

TEST_CASE("matrix get set sub matrix test", "[matrix]") {
  matrix<matrix_storage_cep<double>> m1(4, 4);
  m1.set_value(1, 1, 1.2);
  m1.set_value(2, 2, 3.4);
  m1.set_value(3, 3, 1.56);
  m1.set_value(3, 4, 1.12);
  m1.set_value(4, 4, 5.4);

  SECTION("get sub matrix") {
    auto m2 = m1.get_sub_matrix(3,2,3,2);
    REQUIRE(m2.get_row() == 2);
    REQUIRE(m2.get_column() == 2);
    REQUIRE(value_equal(m2.get_value(1, 2), 1.12));
    REQUIRE(value_equal(m2.get_value(2, 2), 5.4));
    REQUIRE(value_equal(m2.get_value(2, 1), 0.0));
  }

  SECTION("set sub matrix") {
    matrix<matrix_storage_cep<double>> m2(3, 3);
    m2.set_value(1, 1, 2.1);
    m2.set_value(1, 2, 0.56);
    m1.set_value_from_matrix(2, 2, m2);
    REQUIRE(value_equal(m1.get_value(2, 2), 2.1));
    REQUIRE(value_equal(m1.get_value(2, 3), 0.56));
    REQUIRE(value_equal(m1.get_value(4, 4), 0.0));
  }
}

TEST_CASE("matrix vector test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 1);
  m1.set_value(1, 1, 1.4);
  m1.set_value(2, 1, 2.1);
  m1.set_value(3, 1, 4.3);
  double norm = m1.get_vector_second_norm();
  REQUIRE(value_equal(norm * norm, 1.4 * 1.4 + 2.1 * 2.1 + 4.3 * 4.3));
  matrix<matrix_storage_cep<double>> m2(1, 3);
  m2.set_value(1, 1, 1.2);
  m2.set_value(1, 3, 6.5);
  REQUIRE(value_equal(m2.get_vector_inner_product(m1), 1.2 * 1.4 + 4.3 * 6.5));
}

TEST_CASE("matrix transposition test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 2);
  m1.set_value(1, 1, 2.3);
  m1.set_value(2, 1, 4.5);
  auto m2 = tr(m1);
  REQUIRE(m2.get_row() == 2);
  REQUIRE(m2.get_column() == 3);
  REQUIRE(value_equal(m2.get_value(1, 1), 2.3));
  REQUIRE(value_equal(m2.get_value(1, 2), 4.5));
}

TEST_CASE("matrix every element test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(3, 3);
  m1.set_value(2, 2, 1.42);
  using iterator = typename matrix<matrix_storage_cep<double>>::column_iterator;
  m1.every_nozero_element([&](iterator it)->void {
    *it += 2;
  });
  REQUIRE(value_equal(m1.get_value(1, 1), 0.0));
  REQUIRE(value_equal(m1.get_value(2, 2), 1.42 + 2));
}

TEST_CASE("matrix element row transform test","[matrix]") {
  matrix<matrix_storage_cep<double>> m1(4, 4);
  m1.set_value(1, 1, 1.02);
  m1.set_value(1, 4, 2.1);
  m1.set_value(2, 3, 3.1);
  m1.set_value(3, 1, 1.04);
  m1.set_value(4, 4, 1.0);

  SECTION("element row transform swap") {
    m1.element_row_transform_swap(2, 3);
    REQUIRE(value_equal(m1.get_value(2, 1), 1.04));
    REQUIRE(value_equal(m1.get_value(3, 3), 3.1));
    REQUIRE(value_equal(m1.get_value(1, 1), 1.02));
  }

  SECTION("element row transform plus") {
    m1.element_row_transform_plus(2, 3, 2.0);
    REQUIRE(value_equal(m1.get_value(2, 3), 3.1));
    REQUIRE(value_equal(m1.get_value(2, 1), 1.04 * 2));
    REQUIRE(value_equal(m1.get_value(3, 1), 1.04));
  }

  SECTION("element row transform multi") {
    m1.element_row_transform_multi(1, 1.5);
    REQUIRE(value_equal(m1.get_value(1, 1), 1.02 * 1.5));
    REQUIRE(value_equal(m1.get_value(1, 4), 2.1 * 1.5));
  }
}
