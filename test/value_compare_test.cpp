#include "../third_party/catch.hpp"
#include "../include/value_compare.h"

using namespace pnmatrix;

TEST_CASE( "value_compare test", "[util]") {
  REQUIRE(value_equal(0.0, 1e-8) == true);
  REQUIRE(value_equal(-1e-8, 0.0) == true);
  REQUIRE(value_equal(0.0, 0.1) == false);
  REQUIRE(value_equal(1.0,1.0 + 1e-5) == false);
}
