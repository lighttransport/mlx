// SPDX-License-Identifier: MIT
// Copyright 2024 - Present: Light Transport Entertainment, Inc.
//
#if defined(NANOCBLAS_COMPAT)
#include "cblas.h"
#endif

#include "debug-macros.hh"

namespace nanocblas {

bool sgemm(const bool opTransA, const bool opTransB, const int M, const int N,
           const int K, const float alpha, const float *A, const int lda,
           const float *B, const int ldb, const float beta, float *C,
           const int ldc) {
  // Initialize C matrix with beta scaling
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      C[i * ldc + j] *= beta;  // Scale C by beta
    }
  }

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      // TODO(LTE): Use double or Kahan sum for better accuracy?
      float sum = 0.0f;
      for (int l = 0; l < K; ++l) {
        float a_val =
            opTransA ? A[l * lda + i] : A[i * lda + l];  // A[i][l] or A[l][i]
        float b_val =
            opTransB ? B[j * ldb + l] : B[l * ldb + j];  // B[l][j] or B[l][j]
        sum += a_val * b_val;
      }
      C[i * ldc + j] += alpha * sum;  // Scale by alpha and add to C
    }
  }

  return true;
}

}  // namespace nanocblas

#if defined(NANOCBLAS_COMPAT)
void cblas_sgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta, float *C,
                 const int ldc) {
  if (Order != CblasRowMajor) {
    DCOUT("ColMajor is not supported");
    return;
  }

  if (!A) {
    DCOUT("A is null");
    return;
  }

  if (!B) {
    DCOUT("B is null");
    return;
  }

  if (!C) {
    DCOUT("C is null");
    return;
  }

  // transA parameter
  //   'N' : op(A) = A
  //   'T' : op(A) = A'
  //   'C' : op(A) = A'
  // same for transB parameter.

  bool opTransA = false;
  if (TransA == CblasNoTrans) {
    opTransA = false;
  } else if (TransA == CblasTrans) {
    opTransA = true;
  } else if (TransA == CblasConjTrans) {
    opTransA = true;
  } else {
    DCOUT("Invalid transA code: " << std::to_string(transA));
    return;
  }

  bool opTransB = false;
  if (TransB == CblasNoTrans) {
    opTransB = false;
  } else if (TransB == CblasTrans) {
    opTransB = true;
  } else if (TransB == CblasConjTrans) {
    opTransA = true;
  } else {
    DCOUT("Invalid transB code: " << std::to_string(transB));
    return;
  }

  if (!nanocblas::sgemm(opTransA, opTransB, M, N, K, alpha, A, lda, B, ldb,
                        beta, C, ldc)) {
    DCOUT("sgemm failed: ");
    return;
  }
}
#endif
