# MathLib

This is my hobby C++ header-only math library including Analysis, Linear Algebra, Signal Processing.
This library is utilized and tested in my project, embedded base-band DSP TMS320C6748 in on-vehicle digital radio.

## 1. Prerequisites

### 1.1. Essential

* C++14 or above compiler

### 1.2. Recommended

Unit test relies on several build tools.
If you just want to use this library, you do not have to run test case and do not have to setup the following tools.

* CMake
* [Google Test](https://github.com/google/googletest)

  Following pages and shell-scripts will help you install Google Test framework.
  * [Google Testを導入してみた](https://qiita.com/y-vectorfield/items/6238cfd2d9c34aefe364)
  * [Setup GoogleTest on Windows.md](https://gist.github.com/motchy869/22d873415722a1c10bc77d3f761339dc)

## 2. Setup

You just have to copy header files into your workspace and include them.

## 3. Build & Run test

On Linux

```sh
./scripts/doDebugBuild.sh
./scripts/doTest.sh
```

On Windows Power Shell

```bat
./scripts/doDebugBuild.bat
./scripts/doTest.bat
```

You can clean output files by running `doClean.sh` / `doClean.bat`, or simply deleting `./build` directory.

## 4. files

|file|description|
|:---|:---|
|include/common.hpp|common settings<ul><li>switch for bug hunting mode (default is OFF)</li></ul>|
|include/analysis.hpp|Analysis library<ul><li>fast version of `sin`, `cos`, `atan2` with 5-coefficients polynomial approximation</li><li>fast square-root using table lookup (for embedded special processor)</li></ul>|
|include/linAlg.hpp|Linear Algebra library<ul><li>real and complex number type are supported</li><li>matrix and vector scaling, addition, multiplication, transpose, conjugate, shaping, etc.</li><li>vector inner/outer product</li><li>LDL decomposition</li><li>linear equation solver (using Gaussian elimination, LDL decomposition)</li><li>functions to print real/complex vector/matrix on console. `printRealVec`, `printComplexVec`, `printRealMat`, `printComplexMat` are useful in debug.</li></ul>|
|include/sigProc.hpp|Signal Processing library<ul><li>convolution</li><li>exponential weighted 1st-order-IIR filter</li></ul>|
|include/RlsFilter.hpp|Recursive Least Square filter|

## 5. Assumptions

* It is assumed that the entries of any matrix are aligned on memory in row-oriented order (so called C-style array).

## 6. Performance option

`include/common.hpp` has a macro constant `MATH_LIB_INLINE_AGGRESSIVELY`, which is defined as `true` in default.
When this macro is defined as `true`, some relatively-small functions are FORCIBLY in-line expanded.
You can find which functions are expanded by searching `#if MATH_LIB_INLINE_AGGRESSIVELY` in MathLib `include` directory.
You may have to switch this macro to `false` when you have to keep the program size small due to lack of memory space.

## 7. Canary Mode

`include/common.hpp` has a macro constant `MATH_LIB_ENABLE_CANARY_MODE`, which is defined as `false` in default.
When this macro is defined as `true`, some functions in MathLib performs costly parameter validations such as range check, and when encounter invalid parameters, the functions print error message to `std::cerr` and exit with `EXIT_FAILURE`.

When the validity of parameters are guaranteed, parameter validations are just waste of time.
It is recommended to enable Canary Mode in experimental stage, and disable when the application becomes stable.
