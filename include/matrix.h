#pragma once
#include "type.h"
#include "operator_proxy.h"
#include "value_compare.h"
#include "matrix_type_traits.h"
#include <functional>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <memory>

namespace pnmatrix {
template<class Container>
class matrix {
public:
  using container_type = Container;
  using value_type = typename container_type::value_type;
  using row_iterator = typename container_type::row_iterator;
  using column_iterator = typename container_type::row_iterator::column_iterator;
  using const_row_iterator = typename container_type::const_row_iterator;
  using const_column_iterator = typename container_type::const_row_iterator::const_column_iterator;

public:
  matrix() = default;

  matrix(size_type row,size_type column):container_(row, column) {}

  ~matrix() = default;

  matrix(const matrix& other):container_(other.container_) {}

  matrix(matrix&& other):container_(std::move(other.container_)) {}

  template <typename Proxy>
  matrix(const Proxy& pro);

  matrix& operator=(const matrix& other) {
    container_ = other.container_;
    return *this;
  }

  matrix& operator=(matrix&& other) {
    container_ = std::move(other.container_);
    return *this;
  }

  bool operator==(const matrix& other) const {
    if (get_row() != other.get_row() || get_column() != other.get_column()) {
      return false;
    }
    for (auto row = begin(); row != end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        if (!value_equal(*col, other.get_value(col.row_index(), col.column_index()))) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const matrix& other) const {
    return !(*this == other);
  }

  inline size_type get_row() const  {
    return container_.get_row();
  }

  inline size_type get_column() const  {
    return container_.get_column();
  }

  matrix get_sub_matrix(size_type row_begin, size_type r, size_type col_begin, size_type c) const {
    matrix result(r, c);
    for (auto row = begin(); row != end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        if (col.row_index() < row_begin || col.column_index() < col_begin) {
          continue;
        }
        if (col.row_index() > row_begin + r - 1 || col.column_index() > col_begin + c - 1) {
          continue;
        }
        size_type rr = col.row_index() - row_begin + 1;
        size_type rc = col.column_index() - col_begin + 1;
        result.set_value(rr, rc, *col);
      }
    }
    return result;
  }

  void set_value(size_type row, size_type column, const value_type& value) {
    range_check(row, column);
    set_value_withoutcheck(row, column, value);
  }

  void set_value_from_matrix(size_type row_begin, size_type column_begin, const matrix& m) {
    range_check(row_begin + m.get_row() - 1, column_begin + m.get_column() - 1);
    for (auto row = begin(); row != end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        if (col.row_index() < row_begin || col.row_index() > row_begin + m.get_row() - 1) {
          continue;
        }
        if (col.column_index() < column_begin || col.column_index() > column_begin + m.get_column() - 1) {
          continue;
        }
        *col = value_type(0);
      }
    }
    for (auto row = m.begin(); row != m.end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        size_type r = row_begin + col.row_index() - 1;
          size_type c = column_begin + col.column_index() - 1;
          set_value_withoutcheck(r, c, *col);
      }
    }
  }
  
  void set_column(size_type column_begin, const matrix& matrix) {
    column_range_check(column_begin);
    assert(matrix.get_column() == 1);
    set_value_from_matrix(1, column_begin, matrix);
  }

  matrix get_nth_column(size_type n) const {
    return get_sub_matrix(1, get_row(), n, 1);
  }

  void add_value(size_type row, size_type column, const value_type& value) {
    range_check(row, column);
    add_value_withoutcheck(row, column, value);
  }

  value_type get_value(size_type row, size_type column) const {
    range_check(row, column);
    return get_value_withoutcheck(row, column);
  }

  row_iterator begin() {
    return (this->container_).begin();
  }

  row_iterator end() {
    return (this->container_).end();
  }

  const_row_iterator begin() const {
    return (this->container_).begin();
  }

  const_row_iterator end() const {
    return (this->container_).end();
  }

  bool inverse_with_ert(matrix& result) {
    assert(get_row() == get_column());
    if (get_row() == 1 && get_column() == 1) {
      result = matrix(1, 1);
      result.set_value(1, 1, 1.0 / get_value(1, 1));
      return true;
    }
    matrix tmpStorage = *this;
    result = get_identity_matrix(get_row());
    for (size_type i = 1; i <= get_row(); ++i) {
      if (value_equal(get_value_withoutcheck(i, i), value_type(0)) == true) {
        bool error = true;
        for (size_type j = i + 1; j <= get_row(); ++j) {
          if (value_equal(get_value_withoutcheck(j, i), value_type(0)) == false) {
            element_row_transform_swap(i, j);
            result.element_row_transform_swap(i, j);
            error = false;
            break;
          }
        }
        if (error == true)
          return false;
      }

      for (size_type j = i + 1; j <= get_row(); ++j) {
        value_type j_i = get_value_withoutcheck(j, i);
        if (value_equal(j_i, value_type(0)) == false) {
          value_type k = value_type(0) - j_i / get_value(i, i);
          element_row_transform_plus(j, i, k);
          result.element_row_transform_plus(j, i, k);
        }
      }
      value_type k = value_type(1.0) / get_value(i, i);
      element_row_transform_multi(i, k);
      result.element_row_transform_multi(i, k);
    }
    for (size_type i = 2; i <= get_row(); ++i) {
      for (size_type j = i - 1; j >= 1; --j) {
        value_type k = value_type(0) - get_value(j, i);
        element_row_transform_plus(j, i, k);
        result.element_row_transform_plus(j, i, k);
      }
    }
    *this = tmpStorage;
    return true;
  }

  value_type get_vector_second_norm() const {
    assert(get_column() == 1);
    value_type sum = value_type(0);
    for (auto row_iter = begin(); row_iter != end(); ++row_iter) {
      for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
        sum += (*colu_iter) * (*colu_iter);
      }
    }
    return std::sqrt(sum);
  }

  value_type get_vector_inner_product(const matrix& m) const {
    assert(get_row() == 1 && m.get_column() == 1);
    assert(get_column() == m.get_row());
    value_type sum = value_type(0);
    const_row_iterator row_iter = begin();
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      sum += *colu_iter * m.get_value_withoutcheck(colu_iter.column_index(), 1);
    }
    return sum;
  }

  size_type get_nth_row_size(size_type row) const {
    row_range_check(row);
    return container_.get_nth_row_size(row);
  }

  void resize(size_type row, size_type column) {
    if (row == get_row() && column == get_column()) {
      return;
    }
    assert(row >= 1 && column >= 1);
    container_.resize(row, column);
  }

  void every_nozero_element(const std::function<void(const_column_iterator iterator)>& func) const {
    for (auto row_iter = begin(); row_iter != end(); ++row_iter) {
      for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
        if(value_equal(*colu_iter, 0.0) == true) {
          continue;
        }
        func(colu_iter);
      }
    }
  }
  
  void every_nozero_element(const std::function<void(column_iterator iterator)>& func) {
    for (auto row_iter = begin(); row_iter != end(); ++row_iter) {
      for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
        if(value_equal(*colu_iter, 0.0) == true) {
          continue;
        }
        func(colu_iter);
      }
    }
  }

  void element_row_transform_swap(size_type row_i, size_type row_j) {
    row_range_check(row_i);
    row_range_check(row_j);
    container_.element_row_transform_swap(row_i, row_j);
  }

  void element_row_transform_multi(size_type row, value_type k) {
    row_range_check(row);
    container_.element_row_transform_multi(row, k);
  }

  void element_row_transform_plus(size_type row_i, size_type row_j, value_type k) {
    row_range_check(row_i);
    row_range_check(row_j);
    container_.element_row_transform_plus(row_i, row_j, k);
  }

  static matrix get_identity_matrix(size_type row) {
    matrix result(row, row);
    for (size_type i = 1; i <= row; ++i) {
      result.set_value_withoutcheck(i, i, value_type(1));
    }
    return result;
  }

  size_type get_element_count() const {
    return container_.get_element_count();
  }

  void delete_row(size_type row) {
    assert(row > 0);
    container_.delete_row(row);
  }

  void delete_column(size_type column) {
    container_.delete_column(column);
  }

private:
  Container container_;

  explicit matrix(Container&& container):container_(std::move(container)) {}

  inline void range_check(size_type row, size_type column) const  {
    assert(row >= 1 && column >= 1 && row <= get_row() && column <= get_column());
  }

  inline void row_range_check(size_type row) const {
    assert(!(row < 1 || row > get_row()));
  }

  inline void column_range_check(size_type column) const {
    assert(!(column < 1 || column > get_column()));
  }

  inline void set_value_withoutcheck(size_type row, size_type column, const value_type& value) {
    container_.set_value(row, column, value);
  }

  inline void add_value_withoutcheck(size_type row, size_type column, const value_type& value) {
    container_.add_value(row, column, value);
  }

  inline value_type get_value_withoutcheck(size_type row, size_type column) const  {
    return container_.get_value(row, column);
  }
};

template<typename Proxy1, typename Proxy2,
  typename std::enable_if<
    both_real<
             std::disjunction<is_op_type<Proxy1>,is_dense_matrix<Proxy1>>::value,
             std::disjunction<is_op_type<Proxy2>,is_dense_matrix<Proxy2>>::value
    >::value, int>::type = 0>
auto operator+(const Proxy1& m1, const Proxy2& m2)->op_add<Proxy1, Proxy2> {
  assert(m1.get_row() == m2.get_row());
  assert(m1.get_column() == m2.get_column());
  return op_add<Proxy1, Proxy2>(m1, m2, m1.get_row(), m1.get_column());
}

template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
auto operator+(const MatrixType& m1, const MatrixType& m2)->MatrixType {
  assert(m1.get_row() == m2.get_row() && m1.get_column() == m2.get_column());
  MatrixType result(m1.get_row(), m1.get_column());
  for (auto row_iter = m1.begin(); row_iter != m1.end(); ++row_iter) {
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      result.set_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
    }
  }
  for (auto row_iter = m2.begin(); row_iter != m2.end(); ++row_iter) {
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      result.add_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
    }
  }
  return result;
}

template<typename Proxy1, typename Proxy2,
  typename std::enable_if<
    both_real<
      std::disjunction<is_op_type<Proxy1>,is_dense_matrix<Proxy1>>::value,
      std::disjunction<is_op_type<Proxy2>,is_dense_matrix<Proxy2>>::value
    >::value, int>::type = 0>
auto operator-(const Proxy1& m1, const Proxy2& m2)->op_sub<Proxy1, Proxy2> {
  assert(m1.get_row() == m2.get_row());
  assert(m1.get_column() == m2.get_column());
  return op_sub<Proxy1, Proxy2>(m1, m2, m1.get_row(), m1.get_column());
}

template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
auto operator-(const MatrixType& m1, const MatrixType& m2)->MatrixType {
  assert(m1.get_row() == m2.get_row() && m1.get_column() == m2.get_column());
  MatrixType result(m1.get_row(), m1.get_column());
  for (auto row_iter = m1.begin(); row_iter != m1.end(); ++row_iter) {
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      result.set_value(colu_iter.row_index(), colu_iter.column_index(), *colu_iter);
    }
  }
  for (auto row_iter = m2.begin(); row_iter != m2.end(); ++row_iter) {
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      result.add_value(colu_iter.row_index(), colu_iter.column_index(), - *colu_iter);
    }
  }
  return result;
}

template<typename Proxy1, typename Proxy2,
  typename std::enable_if<
    both_real<
      std::disjunction<is_op_type<Proxy1>,is_dense_matrix<Proxy1>>::value,
      std::disjunction<is_op_type<Proxy2>,is_dense_matrix<Proxy2>>::value
    >::value, int>::type = 0>
auto operator*(const Proxy1& m1, const Proxy2& m2)->op_mul<Proxy1, Proxy2> {
  assert(m1.get_column() == m2.get_row());
  return op_mul<Proxy1, Proxy2>(m1, m2, m1.get_row(), m2.get_column());
}

template<typename MatrixType, typename MatrixType2,
  typename std::enable_if<
    std::conjunction<is_sparse_matrix<MatrixType>,is_matrix_type<MatrixType2>>::value,
  int>::type = 0>
auto operator*(const MatrixType& m1, const MatrixType2& m2)->MatrixType2 {
  assert(m1.get_column() == m2.get_row());
  using value_type = typename MatrixType::value_type;
  static_assert (std::is_same<value_type, typename MatrixType2::value_type>::value,"error.");

  MatrixType2 result(m1.get_row(), m2.get_column());
  for (auto row = m1.begin(); row!= m1.end(); ++row) {
    for (size_type i = 1; i <= m2.get_column(); ++i) {
      value_type sum = value_type(0);
      size_type row_ = row.row_index();
      size_type colu_ = i;
      for (auto col = row.begin(); col != row.end(); ++col) {
        sum += *col * m2.get_value(col.column_index(), i);
      }
      result.set_value(row_, colu_, sum);
    }
  }
  return result;
}

template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
MatrixType operator/(const MatrixType& m, typename MatrixType::value_type value) {
  MatrixType result(m.get_row(), m.get_column());
  m.every_nozero_element([&](typename MatrixType::const_column_iterator iter)->void {
    result.set_value(iter.row_index(), iter.column_index(), *iter / value);
  });
  return result;
}

template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
MatrixType operator*(const MatrixType& m, typename MatrixType::value_type value) {
  MatrixType result(m.get_row(), m.get_column());
  m.every_nozero_element([&] (typename MatrixType::const_column_iterator iter)->void {
    result.set_value(iter.row_index(), iter.column_index(), *iter * value);
  });
  return result;
}

template<typename Proxy, typename std::enable_if<
  std::disjunction<is_op_type<Proxy>,is_dense_matrix<Proxy>>::value, int>::type = 0>
auto operator/(const Proxy& m, typename Proxy::value_type value)->op_div_value<Proxy> {
  return op_div_value<Proxy>(m, value, m.get_row(), m.get_column());
}

template<typename Proxy, typename std::enable_if<
  std::disjunction<is_op_type<Proxy>,is_dense_matrix<Proxy>>::value, int>::type = 0>
auto operator*(const Proxy& m, typename Proxy::value_type value)->op_mul_value<Proxy> {
  return op_mul_value<Proxy>(m, value, m.get_row(), m.get_column());
}

template<typename Proxy, typename std::enable_if<
  std::disjunction<is_op_type<Proxy>,is_dense_matrix<Proxy>>::value, int>::type = 0>
op_tr<Proxy> tr(const Proxy& m) {
  return op_tr<Proxy>(m, m.get_column(), m.get_row());
}

template<typename MatrixType, typename std::enable_if<is_sparse_matrix<MatrixType>::value, int>::type = 0>
MatrixType tr(const MatrixType& m) {
  MatrixType result(m.get_column(), m.get_row());
  for (auto row_iter = m.begin(); row_iter != m.end(); ++row_iter) {
    for (auto colu_iter = row_iter.begin(); colu_iter != row_iter.end(); ++colu_iter) {
      result.set_value(colu_iter.column_index(), colu_iter.row_index(), *colu_iter);
    }
  }
  return result;
}

template <typename MatrixType, typename Proxy,typename std::enable_if<
  is_dense_matrix<MatrixType>::value, int>::type = 0>
void construct_from_proxy(MatrixType& self, const Proxy& pro) {
  for(int i = 1; i <= self.get_row(); ++i) {
    for(int j = 1; j <= self.get_column(); ++j) {
      self.set_value(i, j, pro.get_value(i, j));
    }
  }
}

template <typename Container>
template <typename Proxy>
matrix<Container>::matrix(const Proxy& pro):container_(pro.get_row(), pro.get_column()) {
  construct_from_proxy(*this, pro);
}

}
