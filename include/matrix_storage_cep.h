#pragma once
#include "type.h"
#include "value_compare.h"
#include "matrix_type_traits.h"
#include "matrix_storage_cep_config.h"
#include <vector>
#include <cassert>
#include <algorithm>

namespace pnmatrix {
template<class ValueType>
class matrix_storage_cep : public sparse_container {
public:
  using value_type = ValueType;

private:
  using self = matrix_storage_cep;

  struct node_ {
    size_type column_;
    value_type value_;
    node_(size_type c, value_type v) :column_(c), value_(v) {}
  };

  struct each_row_container_ {
      std::vector<node_> this_row_;
  };

public:
  matrix_storage_cep(size_type row, size_type column):
                      my_row_(row),
                      my_column_(column),
                      element_count_(0) {
    assert(row > 0 && column > 0);
    container_.resize(row + 1); // start by 1.

  }
  ~matrix_storage_cep() = default;
  matrix_storage_cep(const self&) = default;
  matrix_storage_cep(self&&) = default;
  self& operator=(const self&) = default;
  self& operator=(self&&) = default;

  bool operator==(const self& other) const {
    if (get_row() != other.get_row() || get_column() != other.get_column()) {
      return false;
    }
    for (auto row = begin(); row != end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        bool e = value_equal(*col, other.get_value(col.row_index(), col.column_index()));
        if (e == false) {
          return false;
        }
      }
    }
    for (auto row = other.begin(); row != other.end(); ++row) {
      for (auto col = row.begin(); col != row.end(); ++col) {
        bool e = value_equal(*col, get_value(col.row_index(), col.column_index()));
        if(e == false) {
          return false;
        }
      }
    }
    return true;
  }

#ifdef DELETE_ZERO
  void set_value(size_type row, size_type column, const value_type& value) {
    std::vector<node_>& row_root = get_nth_row(row);
    bool iszero = value_equal(value, value_type(0));
    bool insert = false;
    bool update = false;
    for (auto it = row_root.begin(); it != row_root.end(); ++it) {
      if (it->column_ == column) {
        it->value_ = value;
        update = true;
        if (iszero) {
          row_root.erase(it);
          --element_count_;
        }
        break;
      }
      else if (it->column_ > column) {
        insert = true;
        break;
      }
    }
    if (update == true || iszero == true) {
      return;
    }
    else if (update == false && insert == false) {
      row_root.push_back(node_(column, value));
      ++element_count_;
      return;
    }
    else {
      row_root.push_back(node_(column, value));
      //排序
      std::sort(row_root.begin(), row_root.end(),
        [](const node_& n1, const node_& n2)->bool {
          return n1.column_ < n2.column_;
      });
      ++element_count_;
      return;
    }
  }
#else
  void set_value(size_type row, size_type column, const value_type& value) {
    std::vector<node_>& row_root = get_nth_row(row);
    auto it = std::find_if(row_root.begin(), row_root.end(), [column](const node_& node)->bool {
      return node.column_ == column;
    });
    if (it == row_root.end()) {
      row_root.push_back(node_(column, value));
      std::sort(row_root.begin(),row_root.end(),[](const node_& n1, const node_& n2)->bool {
        return n1.column_ < n2.column_;
      });
      ++element_count_;
    }
    else {
      it->value_ = value;
    }
  }
#endif

#ifdef DELETE_ZERO
  void add_value(size_type row, size_type column, const value_type& value) {
    bool iszero = value_equal(value, value_type(0));
    if (iszero == true)
      return;
    std::vector<node_>& row_root = get_nth_row(row);
    bool insert = false;
    bool update = false;
    for (auto it = row_root.begin(); it != row_root.end(); ++it) {
      if (it->column_ == column) {
        it->value_ += value;
        update = true;
        if (value_equal(it->value_,value_type(0))) {
          row_root.erase(it);
          --element_count_;
        }
        break;
      }
      else if (it->column_ > column) {
        insert = true;
        break;
      }
    }
    if (update == true) {
      return;
    }
    else if (update == false && insert == false) {
      row_root.push_back(node_(column, value));
      ++element_count_;
      return;
    }
    else {
      row_root.push_back(node_(column, value));
      //排序
      std::sort(row_root.begin(), row_root.end(),
        [](const node_& n1, const node_& n2)->bool {
        return n1.column_ < n2.column_;
      });
      ++element_count_;
      return;
    }
  }
#else
  void add_value(size_type row, size_type column, const value_type& value) {
    std::vector<node_>& row_root = get_nth_row(row);
    auto it = std::find_if(row_root.begin(), row_root.end(), [column](const node_& node)->bool {
      return node.column_ == column;
    });
    if (it == row_root.end()) {
     row_root.push_back(node_(column, value));
     std::sort(row_root.begin(),row_root.end(),[](const node_& n1, const node_& n2)->bool {
       return n1.column_ < n2.column_;
     });
     ++element_count_;
    }
    else {
      it->value_ += value;
    }
  }
#endif

  value_type get_value(size_type row, size_type column) const  {
    for (auto it = get_nth_row(row).begin(); it != get_nth_row(row).end(); ++it) {
      if (it->column_ == column)
        return it->value_;
    }
    return value_type(0);
  }

  inline size_type get_row() const {
    return my_row_;
  }
  
  inline size_type get_column() const   {
    return my_column_;
  }

  size_type get_nth_row_size(size_type row) const {
    return get_nth_row(row).size();
  }

  void delete_row(size_type row) {
    element_count_ -= get_nth_row_size(row);
    auto iter = container_.begin();
    iter += (row);
    container_.erase(iter);
    --my_row_;
  }

  void delete_column(size_type column) {
    --my_column_;
    for (auto each_row = container_.begin(); each_row != container_.end(); ++each_row) {
      auto& this_row = each_row->this_row_;
      for (auto col = this_row.begin(); col != this_row.end();) {
        if (col->column_ == column) {
          col = this_row.erase(col);
          --element_count_;
        }
        else {
          if (col->column_ > column) {
            col->column_ -= 1;
          }
          ++col;
        }
      }
    }
  }

  void resize(size_type new_row, size_type new_column) {
    assert(new_row > 0 && new_column > 0);
    if (new_row != my_row_) {
      if (new_row < my_row_) {
        for(size_type i = new_row + 1; i<=my_row_; ++i) {
          element_count_ -= get_nth_row_size(i);
        }
      }
      container_.resize(new_row + 1);
    }
    if (new_column < my_column_) {
      for (auto row_iter = container_.begin(); row_iter != container_.end(); ++ row_iter) {
        for (auto col_iter = row_iter->this_row_.begin(); col_iter != row_iter->this_row_.end();) {
          if (col_iter->column_ > new_column) {
            col_iter = row_iter->this_row_.erase(col_iter);
            --element_count_;
          }
          else {
            ++col_iter;
          }
        }
      }
    }
    my_row_ = new_row;
    my_column_ = new_column;
  }

  void element_row_transform_swap(size_type row_i, size_type row_j) {
    std::swap(get_nth_row(row_i), get_nth_row(row_j));
  }

  void element_row_transform_multi(size_type row, value_type k) {
    std::vector<node_>& this_row = get_nth_row(row);
    for (auto colu_iter = this_row.begin(); colu_iter != this_row.end(); ++colu_iter) {
      colu_iter->value_ = colu_iter->value_ * k;
    }
  }

  void element_row_transform_plus(size_type row_i, size_type row_j, value_type k) {
    std::vector<node_>& this_row = get_nth_row(row_j);
    for (auto colu_iter = this_row.begin(); colu_iter != this_row.end(); ++colu_iter) {
      size_type column = colu_iter->column_;
      add_value(row_i, column, colu_iter->value_ * k);
    }
  }

  size_type get_element_count() const {
    return element_count_;
  }

  class row_iterator {
  private:
    matrix_storage_cep<value_type>* handle_;
    typename std::vector<each_row_container_>::iterator proxy_;
    size_type row_index_;

  public:
    row_iterator(matrix_storage_cep<ValueType>* h, typename std::vector<each_row_container_>::iterator it, size_type r):
        handle_(h),
        proxy_(it),
        row_index_(r)  {

    }

    row_iterator& operator++() {
      ++row_index_;
      ++proxy_;
      return *this;
    }

    row_iterator operator++(int) {
      row_iterator result = *this;
      ++ *this;
      return result;
    }

    bool operator==(const row_iterator& other) const {
      return proxy_ == other.proxy_;
    }

    bool operator!=(const row_iterator& other) const {
      return proxy_ != other.proxy_;
    }

    size_type row_index() const {
      return row_index_;
    }

    class column_iterator {
    private:
      typename std::vector<node_>::iterator proxy_;
      size_type row_;

    public:
      column_iterator(typename std::vector<node_>::iterator it, size_type r):proxy_(it),row_(r) {

      }

      column_iterator& operator++() {
        ++proxy_;
        return *this;
      }

      column_iterator operator++(int) {
        column_iterator result = *this;
        ++ *this;
        return result;
      }

      bool operator==(const column_iterator& other) const {
        return proxy_ == other.proxy_;
      }

      bool operator!=(const column_iterator& other) const {
        return proxy_ != other.proxy_;
      }

      value_type& operator*() {
        return proxy_->value_;
      }

      value_type* operator->() {
        return &(proxy_->value);
      }

      size_type column_index() const {
        return proxy_->column_;
      }

      size_type row_index() const {
        return row_;
      }
    };

    column_iterator begin() {
      return column_iterator(handle_->container_.at(row_index_).this_row_.begin(), row_index_);
    }

    column_iterator end() {
      return column_iterator(handle_->container_.at(row_index_).this_row_.end(), row_index_);
    }
  };

  class const_row_iterator {
  private:
    const matrix_storage_cep<ValueType>* const handle_;
    typename std::vector<each_row_container_>::const_iterator proxy_;
    size_type row_index_;

  public:
    const_row_iterator(const matrix_storage_cep<ValueType>*const h, typename std::vector<each_row_container_>::const_iterator it, size_type r):
        handle_(h),
        proxy_(it),
        row_index_(r)  {

    }

    const_row_iterator& operator++() {
      ++row_index_;
      ++proxy_;
      return *this;
    }

    const_row_iterator operator++(int) {
      const_row_iterator result = *this;
      ++ *this;
      return result;
    }

    bool operator==(const const_row_iterator& other) const {
      return proxy_ == other.proxy_;
    }

    bool operator!=(const const_row_iterator& other) const {
      return proxy_ != other.proxy_;
    }

    size_type row_index() const {
      return row_index_;
    }

    class const_column_iterator {
    private:
      typename std::vector<node_>::const_iterator proxy_;
      size_type row_;

    public:
      const_column_iterator(typename std::vector<node_>::const_iterator it, size_type r):proxy_(it),row_(r) {

      }

      const_column_iterator& operator++() {
        ++proxy_;
        return *this;
      }

      const_column_iterator operator++(int) {
        const_column_iterator result = *this;
        ++ *this;
        return result;
      }

      bool operator==(const const_column_iterator& other) const {
        return proxy_ == other.proxy_;
      }

      bool operator!=(const const_column_iterator& other) const {
        return proxy_ != other.proxy_;
      }

      const value_type& operator*() {
        return proxy_->value_;
      }

      const value_type* const operator->() {
        return &(proxy_->value);
      }

      size_type column_index() const {
        return proxy_->column_;
      }

      size_type row_index() const {
        return row_;
      }
    };

    const_column_iterator begin() const {
      return const_column_iterator(handle_->container_.at(row_index_).this_row_.begin(), row_index_);
    }

    const_column_iterator end() const {
      return const_column_iterator(handle_->container_.at(row_index_).this_row_.end(), row_index_);
    }
  };

  row_iterator begin() {
    return row_iterator(this,container_.begin() + 1,1);
  }

  row_iterator end() {
    return row_iterator(this,container_.end(), -1);
  }

  const_row_iterator begin() const {
    return const_row_iterator(this, container_.begin() + 1, 1);
  }

  const_row_iterator end() const {
    return const_row_iterator(this, container_.end(), -1);
  }

private:
  size_type my_row_;
  size_type my_column_;
  size_type element_count_;
  std::vector<each_row_container_> container_;

  inline const std::vector<node_>& get_nth_row(size_type row) const {
    return container_.at(row).this_row_;
  }

  inline std::vector<node_>& get_nth_row(size_type row) {
    return container_.at(row).this_row_;
  }
};
}
