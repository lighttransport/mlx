// SPDX-License-Identifier: MIT and Apache 2.0

#if defined(NANOCBLAS_LAPACK)

#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

#include "debug-macros.hh"

namespace nanolapack {

namespace detail {

void inverse(const float *src, int N, float *dst, float *tempData) {
  // Based on MNN
  // https://github.com/alibaba/MNN/blob/master/source/math/Matrix.cpp

  int i, j, k;
  float maxval, temp;
  float *dstData = dst;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      *(dstData + i * N + j) = (i == j) ? 1.0f : 0.0f;
    }
  }

  for (i = 0; i < N; ++i) {
    maxval = *(tempData + i * N + i);
    k = i;
    for (j = i + 1; j < N; ++j) {
      float data1 = *(tempData + j * N + i);
      if (std::fabs(data1) > std::fabs(maxval)) {
        maxval = data1;
        k = j;
      }
    }
    if (k != i) {
      for (j = 0; j < N; ++j) {
        temp = *(tempData + i * N + j);
        *(tempData + i * N + j) = *(tempData + k * N + j);
        *(tempData + k * N + j) = temp;
        temp = *(dstData + i * N + j);
        *(dstData + i * N + j) = *(dstData + k * N + j);
        *(dstData + k * N + j) = temp;
      }
    }
    temp = *(tempData + i * N + i);

    if (std::fabs(temp) < 1.0e-6f) {
      DCOUT("Determinant is zero");
      return;
    }

    for (j = 0; j < N; ++j) {
      *(tempData + i * N + j) = *(tempData + i * N + j) / temp;
      *(dstData + i * N + j) = *(dstData + i * N + j) / temp;
    }

    for (j = 0; j < N; ++j) {
      if (j != i) {
        temp = *(tempData + j * N + i);
        for (k = 0; k < N; ++k) {
          *(tempData + j * N + k) =
              *(tempData + j * N + k) - *(tempData + i * N + k) * temp;
          *(dstData + j * N + k) =
              *(dstData + j * N + k) - *(dstData + i * N + k) * temp;
        }
      }
    }
  }
}

void lu(int n, float *A, int lda, int *ipiv) {
    int i, j, k;
    float temp = 0.0f;
    float alpha = 0.0f;

    for (k = 0; k < n; k++) {
        // Find the pivot element
        alpha = std::fabs(A[k * lda + k]);
        int pivot_index = k;
        for (i = k + 1; i < n; i++) {
            temp = fabs(A[i * lda + k]);
            if (temp > alpha) {
                alpha = temp;
                pivot_index = i;
            }
        }

        // Swap rows k and pivot_index
        ipiv[k] = pivot_index;
        if (pivot_index != k) {
            for (j = 0; j < n; j++) {
                temp = A[k * lda + j];
                A[k * lda + j] = A[pivot_index * lda + j];
                A[pivot_index * lda + j] = temp;
            }
        }

        // Compute the multipliers
        if (std::fabs(A[k * lda + k]) < 1.0e-6f) {
            for (i = k + 1; i < n; i++) {
                A[i * lda + k] /= A[k * lda + k];
            }
        }

        // Update the trailing submatrix
        for (j = k + 1; j < n; j++) {
            for (i = k + 1; i < n; i++) {
                A[i * lda + j] -= A[i * lda + k] * A[k * lda + j];
            }
        }
    }
}

void sgetri(int n, float *A, int lda, int *ipiv) {
    int i, j, k;
    float temp;

    // Initialize the inverse matrix to the identity matrix
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                A[i * lda + j] = 1.0;
            } else {
                A[i * lda + j] = 0.0;
            }
        }
    }

    // Compute the inverse using the LU factorization
    for (k = 0; k < n; k++) {
        // Apply the permutation matrix P
        for (i = 0; i < n; i++) {
            temp = A[i * lda + k];
            A[i * lda + k] = A[ipiv[k] * lda + k];
            A[ipiv[k] * lda + k] = temp;
        }

        // Apply the lower triangular matrix L
        for (i = k + 1; i < n; i++) {
            temp = A[i * lda + k];
            for (j = k; j < i; j++) {
                A[i * lda + j] -= A[i * lda + k] * A[k * lda + j];
            }
            A[i * lda + k] = temp;
        }

        // Apply the upper triangular matrix U
        for (i = k; i < n; i++) {
            temp = A[k * lda + i];
            for (j = k; j < i; j++) {
                A[k * lda + i] -= A[k * lda + j] * A[j * lda + i];
            }
            A[k * lda + i] = temp;
        }
    }

}

// QR factorization
template<typename T>
void qrf(int m, int n, T *A, T *Q, T *R) {
  std::vector<T> v;
  v.resize(m);

    for (int k = 0; k < n && k < m; k++) {
        // Compute the 2-norm of the k-th column
        double norm = 0;
        for (int i = k; i < m; i++) {
            norm += A[i * n + k] * A[i * n + k];
        }
        norm = std::sqrt(norm);

        double alpha = A[k * n + k] > 0 ? -norm : norm;
        double r = std::sqrt(static_cast<T>(0.5) * (norm * norm - A[k * n + k] * alpha));

        // Create the Householder vector
        //double *v = (double *)malloc(m * sizeof(double));
        v[k] = (A[k * n + k] - alpha) / (2 * r);
        for (int i = k + 1; i < m; i++) {
            v[i] = A[i * n + k] / (2 * r);
        }

        // Apply the Householder transformation
        for (int j = k; j < n; j++) {
            double sum = 0;
            for (int i = k; i < m; i++) {
                sum += v[i] * A[i * n + j];
            }
            for (int i = k; i < m; i++) {
                A[i * n + j] -= sum * v[i];
            }
        }

        // Store the results in R
        for (int i = k; i < m && i < n; i++) {
            R[k * n + i] = A[k * n + i];
        }

        // Store Q as the identity matrix with Householder reflections
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j) {
                    Q[i * m + j] = 1.0;
                } else {
                    Q[i * m + j] = 0.0;
                }
            }
        }
        // Update Q with the Householder transformation
        for (int j = 0; j < m; j++) {
            double sum = 0;
            for (int i = k; i < m; i++) {
                sum += v[i] * Q[i * m + j];
            }
            for (int i = k; i < m; i++) {
                Q[i * m + j] -= sum * v[i];
            }
        }

    }
}


}  // namespace detail

void inverse(const float *src, int N, float *dst) {
  std::vector<float> tempMat;
  tempMat.resize(N * N);
  memcpy(tempMat.data(), src, N * N * sizeof(float));

  detail::inverse(src, N, dst, tempMat.data());
}

void qrff(int m, int n, float *A, float *Q, float *R) {
  detail::qrf(m, n, A, Q, R);
}

void qrfd(int m, int n, double *A, double *Q, double *R) {
  detail::qrf(m, n, A, Q, R);
}

#if defined(NANOCBLAS_COMPAT)

void sgetrf_(int n, float *A, int lda, int *ipiv) {
  detail::lu(n, A, lda, ipiv);
}

void sgetri_(int n, float *A, int lda, int *ipiv) {
  detail::sgetri(n, A, lda, ipiv);
}

#endif

}  // namespace nanolapack

#endif // NANOCBAS_LAPACK
