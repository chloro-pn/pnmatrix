#include "../include/matrix.h"
#include "../include/matrix_storage_cep.h"
#include "../include/gmres_solver.h"
#include "../include/gauss_seidel_solver.h"
#include "../include/jacobian_solver.h"
#include <iostream>

using namespace pnmatrix;

template< typename MatrixType>
void print_matrix(const MatrixType& m) {
  for(auto row = m.begin(); row != m.end(); ++row) {
    for(auto col = row.begin(); col != row.end(); ++col) {
      std::cout << "(" << col.row_index() << ", " << col.column_index() << ") " << *col << " ";
    }
    std::cout << "\n";
  }
}

void gmres_example() {
  matrix<matrix_storage_cep<double>> m(3, 3);
  m.set_value(1, 1, 1);
  m.set_value(1, 2, 1);
  m.set_value(1, 3, 1);
  m.set_value(2, 1, 0);
  m.set_value(2, 2, 4);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 2);
  m.set_value(3, 2, -2);
  m.set_value(3, 3, 1);
  std::cout << "example grmes.\n";
  std::cout << "matrix A : \n";
  print_matrix(m);
  matrix<matrix_storage_cep<double>> b(3, 1);
  b.set_value(1, 1, 6);
  b.set_value(2, 1, 5);
  b.set_value(3, 1, 1);

  std::cout << "matrix b : \n";
  print_matrix(b);
  std::cout << "use gmres method to solve A * x = b : \n";
  gmres::option op;
  op.m = 140;
  op.rm = 1e-6;
  std::cout << "restart m : " << op.m << ", error : " << op.rm << std::endl;
  gmres solver(op);
  auto result = solver.solve(m, b);
  std::cout << "result : x \n";
  print_matrix(result);
  std::cout << "###\n";
}

void jacobian_example() {
  matrix<matrix_storage_cep<double>> m(3, 3);
  m.set_value(1, 1, 8);
  m.set_value(1, 2, -3);
  m.set_value(1, 3, 2);
  m.set_value(2, 1, 4);
  m.set_value(2, 2, 11);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 6);
  m.set_value(3, 2, 3);
  m.set_value(3, 3, 12);
  std::cout << "example jacobian.\n";
  std::cout << "matrix A : \n";
  print_matrix(m);
  matrix<matrix_storage_cep<double>> b(3, 1);
  b.set_value(1, 1, 20);
  b.set_value(2, 1, 33);
  b.set_value(3, 1, 36);
  std::cout << "matrix b : \n";
  print_matrix(b);
  std::cout << "use jacobian method to solve A * x = b : \n";
  jacobian::option op;
  op.rm = 1e-6;
  std::cout << "rm : "<< op.rm << std::endl;
  jacobian solver(op);
  std::cout << "result : x \n";
  auto result = solver.solve(m, b);
  print_matrix(result);
  std::cout << "###\n";
}

void gauss_seidel_example() {
  matrix<matrix_storage_cep<double>> m(3, 3);
  m.set_value(1, 1, 8);
  m.set_value(1, 2, -3);
  m.set_value(1, 3, 2);
  m.set_value(2, 1, 4);
  m.set_value(2, 2, 11);
  m.set_value(2, 3, -1);
  m.set_value(3, 1, 6);
  m.set_value(3, 2, 3);
  m.set_value(3, 3, 12);
  std::cout << "example gauss seidel.\n";
  std::cout << "matrix A : \n";
  print_matrix(m);
  matrix<matrix_storage_cep<double>> b(3, 1);
  b.set_value(1, 1, 20);
  b.set_value(2, 1, 33);
  b.set_value(3, 1, 36);
  std::cout << "matrix b : \n";
  print_matrix(b);
  std::cout << "use gauss seidel method to solve A * x = b : \n";
  gauss_seidel::option op;
  op.rm = 1e-6;
  std::cout << "rm : "<< op.rm << std::endl;
  gauss_seidel solver(op);
  std::cout << "result : x \n";
  auto result = solver.solve(m, b);
  print_matrix(result);
  std::cout << "###\n";
}

int main(int argc, char* argv[]) {
  gmres_example();
  jacobian_example();
  gauss_seidel_example();
  return 0;
}
