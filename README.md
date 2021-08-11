# MathLib

My hobby math library includes Linear Algebra, Signal Processing and so on.

## 1. Prerequisites

[Google Test](https://github.com/google/googletest) is mandatory.
Following page and shell-script will help you install Google Test framework.

* [Google Testを導入してみた](https://qiita.com/y-vectorfield/items/6238cfd2d9c34aefe364)
* [Setup GoogleTest on Windows.md](https://gist.github.com/motchy869/22d873415722a1c10bc77d3f761339dc)

## 2. Build & Run test

On Linux

```sh
./doDebugBuild.sh
./doTest.sh
```

On Windows

```bat
./doDebugBuild.bat
./doTest.bat
```

You can clean output files by running `doClean.sh` / `doClean.bat`.

## 3. files

|file|description|
|:---|:---|
|include/analysis.hpp|Analysis library<ul><li>fast version of `sin`, `cos`, `atan2` with 5-coefficients polynomial approximation</li></ul>|
|include/linAlg.hpp|Linear Algebra library<ul><li>real and complex number type are supported</li><li>matrix and vector addition, multiplication</li><li>vector inner product</li><li>LDL decomposition</li><li>linear equation solver (using Gaussian elimination, LDL decomposition)</li></ul>|
|include/sigProc.hpp|Signal Processing library<ul><li>convolution</li><li>exponential weighted 1degree-IIR filter</li></ul>|
|include/RlsFilter.hpp|Recursive Least Square filter|
