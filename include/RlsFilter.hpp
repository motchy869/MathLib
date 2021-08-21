/**
 * @author motchy (motchy869[at]gmail.com)
 * @brief Recursive Least Square filter
 */
#ifndef __MATH_LIB__RLS_FILTER_HPP__
#define __MATH_LIB__RLS_FILTER_HPP__

#define ENVIRONMENT_STANDARD 0
#define ENVIRONMENT_TI_RTOS  1
#define ENVIRONMENT ENVIRONMENT_STANDARD

#include <cassert>
#include <cstddef>
#include <exception>
#include <new>
#include <stdexcept>
#if ENVIRONMENT == ENVIRONMENT_TI_RTOS
    #include <xdc/runtime/Error.h>
    #include <xdc/runtime/IHeap.h>
    #include <xdc/runtime/Memory.h>
#endif
#include "common.hpp"
#include "linAlg.hpp"
#include "sigProc.hpp"

namespace MathLib {
    /**
     * @brief An extended implementation of RLS filter.
     * @details see: https://en.wikipedia.org/wiki/Recursive_least_squares_filter @n
     * The extensions in this code are: @n
     *   - complex number data
     *   - "thinned" (or "decimated")  estimation
     *
     * @tparam T the number type of RLS filter input.
     */
    template <typename T>
    class RlsFilter {
        private:
            const int m_p;
            const int m_thinningRate;

            /* working space */
            char *m_workspace; // used for `m_R_x`, ..., `m_E_solveLinEq`
            size_t m_bufSize_workspace; // the capacity in bytes of `m_workspace`
            T *m_R_x;
            T *m_r_d;
            T *m_x;
            T *m_x_conj;
            T *m_delta_R_x;
            T *m_delta_r_d;
            char *m_workspace_solveLinEq; // workspace for `MathLib::LinAlg::solveLinearEquation` function

        public:
            /**
             * @brief Constructs a new Rls Filter object
             *
             * @param[in] p the filter order, must be 1 or greater.
             * @param[in] thinningRate thinningRate The input/output data rate, defaults to 1. When `thinningRate` >= 2, it is assumed that input is over-sampled at rate `thinningRate`. This value must be in range [1, p].
             */
            RlsFilter(int p, int thinningRate=1): m_p(p), m_thinningRate(thinningRate) {
                #if ENVIRONMENT == ENVIRONMENT_STANDARD
                    if (p < 1) {
                        throw std::out_of_range("`p` must be 1 or greater.");
                    }
                    if (thinningRate < 1) {
                        throw std::out_of_range("`thinningRate` must be 1 or greater.");
                    }
                    if (p < thinningRate) {
                        throw std::out_of_range("`p` must NOT be less than `thinningRate`.");
                    }
                #elif ENVIRONMENT == ENVIRONMENT_TI_RTOS
                        if (p < 1) {
                            System_abort("`p` must be 1 or greater.");
                        }
                        if (thinningRate < 1) {
                            System_abort("`thinningRate` must be 1 or greater.");
                        }
                        if (p < thinningRate) {
                            System_abort("`p` must NOT be less than `thinningRate`.");
                        }
                #endif

                /* Allocates memory for working space. */
                const size_t bufSize_R_x = (p+1)*(p+1)*sizeof(T);
                const size_t bufSize_r_d = (p+1)*sizeof(T);
                const size_t bufSize_x = (p+1)*sizeof(T);
                const size_t bufSize_x_conj = bufSize_x;
                const size_t bufSize_delta_R_x = bufSize_R_x;
                const size_t bufSize_delta_r_d = bufSize_r_d;
                const size_t bufSize_workspace_solveLinEq = (p+1)*(p+2)*sizeof(T) + (p+1)*sizeof(size_t);
                m_bufSize_workspace = bufSize_R_x + bufSize_r_d + bufSize_x + bufSize_x_conj + bufSize_delta_R_x + bufSize_delta_r_d + bufSize_workspace_solveLinEq;

                #if ENVIRONMENT == ENVIRONMENT_STANDARD
                    m_workspace = reinterpret_cast<char *>(malloc(m_bufSize_workspace));
                #elif ENVIRONMENT == ENVIRONMENT_TI_RTOS
                    Error_Block eb;
                    Error_init(&eb);
                    m_workspace = reinterpret_cast<char *>(Memory_alloc(NULL, m_bufSize_workspace, 0, &eb));
                #endif
                if (nullptr == m_workspace) {
                    #if ENVIRONMENT == ENVIRONMENT_STANDARD
                        throw std::bad_alloc();
                    #elif ENVIRONMENT == ENVIRONMENT_TI_RTOS
                        System_abort("failed to allocate memory.");
                    #endif
                }

                char *ptr = m_workspace;
                m_R_x                  = reinterpret_cast<T *>(ptr); ptr += bufSize_R_x;
                m_r_d                  = reinterpret_cast<T *>(ptr); ptr += bufSize_r_d;
                m_x                    = reinterpret_cast<T *>(ptr); ptr += bufSize_x;
                m_x_conj               = reinterpret_cast<T *>(ptr); ptr += bufSize_x_conj;
                m_delta_R_x            = reinterpret_cast<T *>(ptr); ptr += bufSize_delta_R_x;
                m_delta_r_d            = reinterpret_cast<T *>(ptr); ptr += bufSize_delta_r_d;
                m_workspace_solveLinEq = ptr;                        ptr += bufSize_workspace_solveLinEq;
                assert(ptr == m_workspace + m_bufSize_workspace);
            }

            ~RlsFilter() {
                if (nullptr != m_workspace) {
                    #if ENVIRONMENT == ENVIRONMENT_STANDARD
                        free(m_workspace);
                    #elif ENVIRONMENT == ENVIRONMENT_TI_RTOS
                        Memory_free(NULL, m_workspace, m_bufSize_workspace);
                    #endif
                }
            }

            /**
             * @brief Calculates an optimal weight vector which estimates `vec_d` from `vec_x`.
             *
             * @param[in] vec_d reference data
             * @param[in] vec_x distorted data
             * @param[in] Ld the length of `vec_d`, must be 1 or greater
             * @param[out] w_opt array of length p+1, the place to put optimal filter coefficient
             *
             * @attention vec_x[-p], ..., vec_x[-1] must be valid data.
             * @attention When `Ld < 1`, this function does nothing and returns immediately in release build, aborts in debug build.
             */
            void trainHard(const T *vec_d, const T *vec_x, int Ld, T *w_opt) {
                #define IS_DEBUGGING false
                assert(Ld >= 1);
                if (Ld < 1) {
                    return;
                }

                constexpr T ZERO = static_cast<T>(0);
                const int r_T = m_thinningRate;
                std::fill_n(m_R_x, (m_p+1)*(m_p+1), ZERO);
                std::fill_n(m_r_d, m_p+1, ZERO);

                #if IS_DEBUGGING == true
                    const char *format = "%9.6f + %9.6fj";
                    puts("before loop:");
                    puts("m_R_x:");
                    MathLib::LinAlg::Debug::printComplexMat(m_p+1, m_p+1, m_R_x, format);
                    puts("m_r_d:");
                    MathLib::LinAlg::Debug::printComplexVec(m_p+1, m_r_d, format);
                    puts("-----");
                #endif

                for (int n=0; n<Ld; ++n) {
                    using namespace MathLib::LinAlg;
                    std::reverse_copy(&vec_x[n*r_T - (Ld-1)], &vec_x[n*r_T+1], m_x);
                    conjugateVec(m_p+1, m_x, m_x_conj);
                    mulMat(m_p+1, 1, m_p+1, m_x_conj, m_x, m_delta_R_x);
                    addMat_inplace(m_p+1, m_p+1, m_R_x, m_delta_R_x);
                    scaleVec(m_p+1, vec_d[n], m_x_conj, m_delta_r_d);
                    addVec_inplace(m_p+1, m_r_d, m_delta_r_d);

                    #if IS_DEBUGGING == true
                        printf("n == %d:\n", n);
                        puts("m_x:");
                        MathLib::LinAlg::Debug::printComplexVec(m_p+1, m_x, format);
                        puts("-----");
                    #endif
                }
                MathLib::LinAlg::solveLinearEquation(m_p+1, m_R_x, m_r_d, w_opt, m_workspace_solveLinEq);

                #if IS_DEBUGGING == true
                    puts("w_opt:");
                    MathLib::LinAlg::Debug::printComplexVec(m_p+1, w_opt, format);
                #endif
                #undef IS_DEBUGGING
            }
    };
}

#undef ENVIRONMENT_STANDARD
#undef ENVIRONMENT_TI_RTOS
#undef ENVIRONMENT

#endif // __MATH_LIB__RLS_FILTER_HPP__
