#pragma once
#include <type_traits>

namespace pnmatrix {
template <bool b1, bool b2>
struct both_real : std::false_type {};

template <>
struct both_real<true, true> : std::true_type {};

class dense_container {
public:
  using dense_tag = void;
};

class sparse_container {
public:
  using sparse_tag = void;
};

template <typename , typename = std::void_t<>>
struct is_op_type : std::false_type {};

template <typename T>
struct is_op_type<T, std::void_t<typename T::op_type_flag>> : std::true_type {};

template <typename T, typename = std::void_t<>>
struct is_dense_matrix : std::false_type {};

template <typename T, typename = std::void_t<>>
struct is_sparse_matrix : std::false_type {};

template <typename T>
class matrix;

template <typename T>
struct is_dense_matrix<matrix<T>, std::void_t<typename T::dense_tag>> : public std::true_type {};

template <typename T>
struct is_sparse_matrix<matrix<T>, std::void_t<typename T::sparse_tag>> : public std::true_type {};

template <typename T>
struct is_matrix_type : std::false_type {};

template <typename T>
struct is_matrix_type<matrix<T>> : std::true_type {};
}
