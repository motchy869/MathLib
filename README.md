# MathLib

My hobby math library including Linear Algebra, Signal Processing and so on.

## 1. Prerequisites

[Google Test](https://github.com/google/googletest) is mandatory.
Following page and shell-script will help you install Google Test framework.

* [Google Testを導入してみた](https://qiita.com/y-vectorfield/items/6238cfd2d9c34aefe364)
* [Setup GoogleTest on Windows.md](https://gist.github.com/motchy869/22d873415722a1c10bc77d3f761339dc)

## 2. Build & Run test

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

You can clean output files by running `doClean.sh` / `doClean.bat`.

## 3. files

|file|description|
|:---|:---|
|include/common.hpp|common settings<ul><li>switch for bug hunting mode (default is OFF)</li></ul>|
|include/analysis.hpp|Analysis library<ul><li>fast version of `sin`, `cos`, `atan2` with 5-coefficients polynomial approximation</li></ul>|
|include/linAlg.hpp|Linear Algebra library<ul><li>real and complex number type are supported</li><li>matrix and vector addition, multiplication</li><li>vector inner/outer product</li><li>LDL decomposition</li><li>linear equation solver (using Gaussian elimination, LDL decomposition)</li><li>functions to print real/complex vector/matrix on console. `printRealVec`, `printComplexVec`, `printRealMat`, `printComplexMat` are useful in debug.</li></ul>|
|include/sigProc.hpp|Signal Processing library<ul><li>convolution</li><li>exponential weighted 1degree-IIR filter</li></ul>|
|include/RlsFilter.hpp|Recursive Least Square filter|

## 4. Assumptions

* It is assumed that the entries of any matrix are aligned on memory in row-oriented order (so called C-style array).

## 5. Bug Hunting Mode

`include/common.hpp` has a macro constant `MATH_LIB_ENABLE_BUG_HUNTING_MODE`, which is defined as `false` in default.
When this macro is defined as `true`, the functions in MathLib performs costly parameter validations such as range check, and when encounter invalid parameters, the functions print error message to `std::cerr` and exit with `EXIT_FAILURE`.

When the validity of parameters are guaranteed, parameter validations are just waste of time.
It is recommended to enable bug hunting mode in experimental stage, and disable when the application becomes stable.
