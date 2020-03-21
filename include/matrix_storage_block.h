#pragma once
#include "matrix_type_traits.h"
#include "value_compare.h"
#include "type.h"
#include <vector>
#include <cassert>

namespace pnmatrix {
template <typename ValueType>
class matrix_storage_block : public dense_container {
public:
  using value_type = ValueType;

  matrix_storage_block(size_type row, size_type column):my_row_(row), my_column_(column),block_row_(row),block_column_(column) {
    assert(row > 0 && column > 0);
    block_.resize(row * column, ValueType(0));
  }

  matrix_storage_block(const matrix_storage_block&) = default;
  matrix_storage_block(matrix_storage_block&&) = default;
  matrix_storage_block& operator=(const matrix_storage_block&) = default;
  matrix_storage_block& operator=(matrix_storage_block&&) = default;

  bool operator==(const matrix_storage_block& other) const {
    if(get_row() != other.get_row() || get_column() != other.get_column()) {
      return false;
    }
    for (int i = 1; i <= my_row_; ++i) {
      for (int j = 1; j <= my_column_; ++j) {
        if (value_equal(get_value(i, j), other.get_value(i, j)) == false) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const matrix_storage_block& other) const {
    return !(*this == other);
  }

  void set_value(size_type row, size_type column, const value_type& v) {
    block_[get_index(row, column)] = v;
  }

  void add_value(size_type row, size_type column, const value_type& v) {
    block_[get_index(row, column)] += v;
  }

  value_type get_value(size_type row, size_type column) const {
    return block_[get_index(row, column)];
  }

  inline size_type get_row() const {
    return my_row_;
  }

  inline size_type get_column() const   {
    return my_column_;
  }

  size_type get_nth_row_size(size_type row) const {
    return my_column_;
  }

  void delete_row(size_type row) {
    for (int i = row + 1; i <= my_row_; ++i ) {
      for (int j = 1; j <= my_column_; ++j) {
        set_value(i-1, j, get_value(i,j));
      }
    }
    --my_row_;
  }

  void delete_column(size_type column) {
    for (int j = column + 1; j <= my_column_; ++j) {
      for (int i = 1; i <= my_row_; ++i) {
        set_value(i, j-1, get_value(i,j));
      }
    }
    --my_column_;
  }

  void resize(size_type new_row, size_type new_column) {
    assert(new_row > 0 && new_column > 0);
    matrix_storage_block tmp(new_row, new_column);
    for (int i = 1; i <= my_row_; ++i) {
      for (int j = 1; j <= my_column_; ++j) {
        if(i <= new_row && j <= new_column) {
          tmp.set_value(i, j, get_value(i, j));
        }
      }
    }
    block_.swap(tmp.block_);
    my_row_ = new_row;
    my_column_ = new_column;
    block_row_ = new_row;
    block_column_ = new_column;
  }

  class row_iterator {
  private:
    matrix_storage_block<value_type>* handle_;
    size_type row_index_;

  public:
    row_iterator(matrix_storage_block<ValueType>* h, size_type r):
        handle_(h),
        row_index_(r)  {

    }

    row_iterator& operator++() {
      ++row_index_;
      return *this;
    }

    row_iterator operator++(int) {
      row_iterator result = *this;
      ++ *this;
      return result;
    }

    bool operator==(const row_iterator& other) const {
      return handle_ == other.handle_ && row_index_ == other.row_index_;
    }

    bool operator!=(const row_iterator& other) const {
      return handle_ != other.handle_ || row_index_ != other.row_index_;
    }

    size_type row_index() const {
      return row_index_;
    }

    class column_iterator {
    private:
      matrix_storage_block<value_type>* handle_;
      size_type row_;
      size_type column_;

    public:
      column_iterator(matrix_storage_block<value_type>* h, size_type r, size_type c):handle_(h),row_(r), column_(c) {

      }

      column_iterator& operator++() {
        ++column_;
        return *this;
      }

      column_iterator operator++(int) {
        column_iterator result = *this;
        ++ *this;
        return result;
      }

      bool operator==(const column_iterator& other) const {
        return handle_ == other.handle_ && row_ == other.row_ && column_ == other.column_;
      }

      bool operator!=(const column_iterator& other) const {
        return handle_ != other.handle_ || row_ != other.row_ || column_ != other.column_;
      }

      value_type& operator*() {
        return handle_->block_[handle_->get_index(row_, column_)];
      }

      value_type* operator->() {
        return &(operator*());
      }

      size_type column_index() const {
        return column_;
      }

      size_type row_index() const {
        return row_;
      }
    };

    column_iterator begin() {
      return column_iterator(handle_, row_index_, 1);
    }

    column_iterator end() {
      return column_iterator(handle_, row_index_, handle_->get_column() + 1);
    }
  };

  class const_row_iterator {
  private:
    const matrix_storage_block<ValueType>* const handle_;
    size_type row_index_;

  public:
    const_row_iterator(const matrix_storage_block<ValueType>*const h, size_type r):
        handle_(h),
        row_index_(r)  {

    }

    const_row_iterator& operator++() {
      ++row_index_;
      return *this;
    }

    const_row_iterator operator++(int) {
      const_row_iterator result = *this;
      ++ *this;
      return result;
    }

    bool operator==(const const_row_iterator& other) const {
      return handle_ == other.handle_ && row_index_ == other.row_index_;
    }

    bool operator!=(const const_row_iterator& other) const {
      return handle_ != other.handle_ || row_index_ != other.row_index_;
    }

    size_type row_index() const {
      return row_index_;
    }

    class const_column_iterator {
    private:
      const matrix_storage_block<value_type>* handle_;
      size_type row_;
      size_type column_;

    public:
      const_column_iterator(const matrix_storage_block<value_type>* h, size_type r, size_type c):handle_(h),row_(r),column_(c) {

      }

      const_column_iterator& operator++() {
        ++column_;
        return *this;
      }

      const_column_iterator operator++(int) {
        const_column_iterator result = *this;
        ++ *this;
        return result;
      }

      bool operator==(const const_column_iterator& other) const {
        return handle_ == other.handle_ && row_ == other.row_ && column_ == other.column_;
      }

      bool operator!=(const const_column_iterator& other) const {
        return handle_ != other.handle_ || row_ != other.row_ || column_ != other.column_;
      }

      const value_type& operator*() {
        return handle_->block_[handle_->get_index(row_, column_)];
      }

      const value_type* const operator->() {
        return &(operator*());
      }

      size_type column_index() const {
        return column_;
      }

      size_type row_index() const {
        return row_;
      }
    };

    const_column_iterator begin() const {
      return const_column_iterator(handle_, row_index_, 1);
    }

    const_column_iterator end() const {
      return const_column_iterator(handle_, row_index_, handle_->get_column() + 1);
    }
  };

  row_iterator begin() {
    return row_iterator(this, 1);
  }

  row_iterator end() {
    return row_iterator(this, my_row_ + 1);
  }

  const_row_iterator begin() const {
    return const_row_iterator(this, 1);
  }

  const_row_iterator end() const {
    return const_row_iterator(this, my_row_ + 1);
  }

  void element_row_transform_swap(size_type row_i, size_type row_j) {
    for (int j = 1; j <= my_column_; ++j) {
      std::swap(block_[get_index(row_i, j)], block_[get_index(row_j, j)]);
    }
  }

  void element_row_transform_multi(size_type row, value_type k) {
    for (int j = 1; j <= my_column_; ++j) {
      block_[get_index(row, j)] *= k;
    }
  }

  void element_row_transform_plus(size_type row_i, size_type row_j, value_type k) {
    for (int j = 1; j <= my_column_; ++j) {
      add_value(row_i, j, get_value(row_j, j) * k);
    }
  }

private:
  size_type my_row_;
  size_type my_column_;
  std::vector<ValueType> block_;
  //这两个值的作用是删除行或者列的时候保留block的映射关系。
  size_type block_row_;
  size_type block_column_;

  inline
  size_type get_index(size_type row, size_type column) const {
    return (row - 1) * block_column_ + column - 1;
  }
};

}
