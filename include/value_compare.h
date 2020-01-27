#pragma once

namespace pnmatrix {
template< typename T>
bool value_equal(const T& v1, const T& v2) {
  return (v1 - v2) > -1e-6 && (v1 - v2) < 1e-6;
}
}
