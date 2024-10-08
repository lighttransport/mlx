#pragma once

namespace nanocblas {

// row-major only
bool sgemm(const bool TransA, const bool TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);

} // namespace nanocblas
