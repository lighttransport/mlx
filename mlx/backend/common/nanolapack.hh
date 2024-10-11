#pragma once

#if defined(NANOCBLAS_LAPACK)

namespace nanolapack {

// Compute the inverse of square matrix N x N
void inverse(const float *src, int N, float *dst);

// QR factorization
void qrff(int m, int n, float *A, float *Q, float *R);
void qrfd(int m, int n, double *A, double *Q, double *R);

} // namespace nanolapack

#if defined(NANOCBLAS_COMPAT)

void sgetrf_(int n, float *A, int lda, int *ipiv);
void sgetri_(int n, float *A, int lda, int *ipiv);
void sgeqrf_(int m, int n, float *A, int lda, float *TAU, float *WORK, int LWORK, int *info);
#endif

#endif
