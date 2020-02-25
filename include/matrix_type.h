#pragma once
#include <type_traits>

namespace pnmatrix {
class dense_container {
public:
  using dense_tag = void;
};

class sparse_container {
public:
  using sparse_tag = void;
};
}
