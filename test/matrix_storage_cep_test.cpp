#include "../third_party/catch.hpp"
#include "../include/matrix_storage_cep.h"
#include "../include/value_compare.h"
#include <vector>

using namespace pnmatrix;

TEST_CASE("matrix storage cep set and get value test ","[matrix_container]") {
  matrix_storage_cep<double> m1(3, 3);
  REQUIRE(m1.get_row() == 3);
  REQUIRE(m1.get_column() == 3);
  REQUIRE(value_equal(m1.get_value(1, 2), 0.0));
  m1.set_value(1, 2, 2.45);
  REQUIRE(value_equal(m1.get_value(1, 2), 2.45));
  m1.add_value(1, 2, 1.23);
  REQUIRE(value_equal(m1.get_value(1, 2), 2.45 + 1.23));
  m1.set_value(1, 2, 0.0);
  REQUIRE(value_equal(m1.get_value(1, 2), 0.0));
}

TEST_CASE("matrix storage cep copy and move test", "[matrix_container]") {
  matrix_storage_cep<double> m1(3, 3);
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

TEST_CASE("matrix storage cep delete row and column test","[matrix_container]") {
  matrix_storage_cep<double> m1(3, 3);
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

TEST_CASE("matrix storage cep resize test","[matrix_container]") {
  matrix_storage_cep<double> m1(4, 4);
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

TEST_CASE("matrix storage cep iterator test", "[matrix_container]") {
  matrix_storage_cep<double> m1(4, 4);
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

TEST_CASE("matrix storage cep element row transform test","[matrix_container]") {
  matrix_storage_cep<double> m1(4, 4);
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
