# pnmatrix
a library for solving linear equations based on c++11.

support:

* Gmres-m.

* Jacobi iteration.

* Gauss-seidel iteration.

# license
Use of this code is governed by a MIT license that can be found in the License file.

# build
pnmatrix is a header-only library, so you can just copy the include folder to your project or add include to your project's include_path.

# test and example
pnmatrix uses Catch2(v2.11.1) for unit test, which you can find in : https://github.com/catchorg/Catch2.
You can use CMake to build test and example executables:
```
mkdir build
cd build
cmake ..
make
```
Then you can find executables in folder build/bin. Try to execute them:)
```
cd bin
./pnmatrix_example
./pnmatrix_test
```
you can find source code in folder examples and test
