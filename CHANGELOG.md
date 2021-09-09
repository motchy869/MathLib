# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Added `MATH_LIB_INLINE_AGGRESSIVELY` option.
- Speeded up following functions:
  - `Analysis::`
    - `conj`, `sqAbs`
  - `LinAlg::`
    - `innerProd`, `vecSelfOuterProd`, `mulMat`, `ldlDecomp`, `solveLinearEquation`, `solveLinEqHermitian`, `hermitianInnerProduct`
  - `SigProc::`
    - `convolve`, `convolve_type2`
- Improved include-guard: replaced old-style guard using `#ifndef` with `#pragma once`.
- Added `LinAlg::fillLowTri` function.
- Added `Analysis::atan_polyApprox_deg7` function.
- Renamed `Analysis::atan_polyApprox` to `Analysis::atan_polyApprox_deg9`.
- Renamed `Analysis::atan2_polyApprox` to `Analysis::atan2_polyApprox_deg9`.

## [0.4.0] - 2021/8/27

- Renamed project: `MotchyMathLib` -> `MathLib`.
- Optimized `atan2_polyApprox` function to speed up on Texas Instruments DSP TMS320C6748.
  - There is no effect on modern PC processors because the compilers for them are clever enough to perform above easy hand-optimization.
- Added `vecSelfOuterProd` function.
- Added bug hunting mode.
- Added `dropSubMat` function.
- Added `setDiag` function.
- Added `addDiag` function.
- Speeded up `ldlDecomp` function, and consequently `solveLinEqHermitian` function is speeded up.
- Add in-place version of `conjugateMat`.
- Add `innerProd` function
- Rename `addMat_inplace`, `addVec_inplace` to `addMat`, `addVec`

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
