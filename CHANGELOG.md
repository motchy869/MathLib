# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Added `dropSubMat` function.
- Added bug hunting mode.
- Added `vecSelfOuterProd` function.
- Optimized `atan2_polyApprox` function to speed up on Texas Instruments DSP TMS320C6748.
  - There is no effect on modern PC processors because the compilers for them are clever enough to perform above easy hand-optimization.
- Renamed project: `MotchyMathLib` -> `MathLib`.

## [0.3.0] - 2021/8/10

- Add `solveLinEqHermitian` function
- Add `ldlDecomp` function
- Add Windows build support
- Add cleaning shell script `doClean.sh`

## [0.2.0] - 2021/7/29

- Refine CMakeLists.txt
- Refine include guard macro name for `analysis.hpp`
- Add `ExpWeightedIirFilter` class
- Add `printRealMat()` and `printRealVec()` functions
- Add `sin_polyApprox()` function
- Add `cos_polyApprox()` function

## [0.1.0] - 2021/7/22

- Add `atan2_polyApprox` function
- Refine namespace
