# MotchyMathLib

My hobby math library includes Linear Algebra, Signal Processing and so on.

## 1. Prerequisites

[Google Test](https://github.com/google/googletest) is mandatory.
Following page and shell-script will help you install Google Test framework.

* [Google Testを導入してみた](https://qiita.com/y-vectorfield/items/6238cfd2d9c34aefe364)
* [install-googletest.sh](https://gist.github.com/motchy869/22d873415722a1c10bc77d3f761339dc)

## 2. Build & Run test

```sh
./doDebugBuild.sh
./doTest.sh
```

## 3. files

|file|description|
|:---|:---|
|include/analysis.hpp|Analysis library<ul><li>fast version of `sin`, `cos`, `atan2` with 5-coefficients polynomial approximation</li></ul>|
|include/linAlg.hpp|Linear Algebra library<ul><li>real and complex number type are supported</li><li>matrix and vector addition, multiplication</li><li>vector inner product</li><li>linear equation solver (Gaussian elimination)</li></ul>|
|include/sigProc.hpp|Signal Processing library<ul><li>convolution</li><li>exponential weighted 1degree-IIR filter</li></ul>|
|include/RlsFilter.hpp|Recursive Least Square filter|
