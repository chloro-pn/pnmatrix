#include "type.h"
#include <type_traits>
#include <cassert>

namespace pnmatrix {
class op_base {
private:
  size_type row_;
  size_type col_;

public:
  using op_type_flag = void;
  op_base(size_type row, size_type col):row_(row), col_(col) {

  }

  op_base(const op_base&) = default;

  size_type get_row() const {
    return row_;
  }
  size_type get_column() const {
    return col_;
  }
};

template<typename Left, typename Right>
class op_add : public op_base {
private:
  static_assert (std::is_same<typename Left::value_type, typename Right::value_type>::value,
  "invalid op_add !");
  const Left& l_;
  const Right& r_;

public:
  using value_type = typename Left::value_type;

  op_add(const Left& l, const Right& r, size_type row, size_type col):op_base(row, col), l_(l), r_(r) {

  }

  value_type get_value(size_type row, size_type col) const {
    return l_.get_value(row, col) + r_.get_value(row, col);
  }
};

template<typename Left, typename Right>
class op_mul : public op_base {
private:
  static_assert (std::is_same<typename Left::value_type, typename Right::value_type>::value,
  "invalid op_mul !");
  const Left& l_;
  const Right& r_;

public:
  using value_type = typename Left::value_type;

  op_mul(const Left& l, const Right& r, size_type row, size_type col):op_base(row, col), l_(l), r_(r) {
    assert(l_.get_column() == r_.get_row());
  }

  value_type get_value(size_type row, size_type col) const {
    value_type sum = value_type(0);
    for(size_type i = 1; i <= l_.get_column(); ++i) {
      sum += l_.get_value(row, i) * r_.get_value(i, col);
    }
    return sum;
  }
};

template<typename Left, typename Right>
class op_sub : public op_base {
private:
  static_assert (std::is_same<typename Left::value_type, typename Right::value_type>::value,
  "invalid op_sub !");
  const Left& l_;
  const Right& r_;

public:
  using value_type = typename Left::value_type;

  op_sub(const Left& l, const Right& r, size_type row, size_type col):op_base(row, col), l_(l), r_(r) {

  }

  value_type get_value(size_type row, size_type col) const {
    return l_.get_value(row, col) - r_.get_value(row, col);
  }
};

template<typename Proxy>
class op_div_value : public op_base {
private:
  const Proxy& l_;
  typename Proxy::value_type v_;

public:
  using value_type = typename Proxy::value_type;

  op_div_value(const Proxy& l,const value_type&v, size_type row, size_type col):op_base(row, col), l_(l), v_(v) {

  }

  value_type get_value(size_type row, size_type col) const {
    return l_.get_value(row, col) / v_;
  }
};

template<typename Proxy>
class op_mul_value : public op_base {
private:
  const Proxy& l_;
  typename Proxy::value_type v_;

public:
  using value_type = typename Proxy::value_type;

  op_mul_value(const Proxy& l,const value_type& v, size_type row, size_type col):op_base(row, col), l_(l), v_(v) {

  }

  value_type get_value(size_type row, size_type col) const {
    return l_.get_value(row, col) * v_;
  }
};

template<typename Proxy>
class op_tr : public op_base {
private:
  const Proxy& l_;

public:
  using value_type = typename Proxy::value_type;

  op_tr(const Proxy& l, size_type row, size_type col):op_base(row, col), l_(l) {

  }

  value_type get_value(size_type row, size_type col) const {
    return l_.get_value(col, row);
  }
};
}
