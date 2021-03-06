/*
 * Copyright 2020 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TINY_MATRIXXXX_H
#define TINY_MATRIXXXX_H

#include <stdio.h>

// #include "tiny_spatial_motion_vector.h"
#include "tiny_vector_x.h"
#include "tiny_vector3.h"

/**
 * Represents a matrix with arbitrary number of columns and custom column type.
 */
template <typename TinyScalar, typename TinyConstants,
          template <typename, typename> typename ColumnType>
class TinyMatrixXxX_ {
  typedef ColumnType<TinyScalar, TinyConstants> ColumnVector;

  // columns are stored as vectors
  ColumnVector* m_columns{nullptr};

  inline void allocate() {
    m_columns = new ColumnVector[m_cols];
    for (int i = 0; i < m_cols; ++i) m_columns[i] = ColumnVector(m_rows);
    set_zero();
  }

 public:
  int m_rows{0};
  int m_cols{0};

  TinyMatrixXxX_() = default;

  TinyMatrixXxX_(int rows, int cols) : m_rows(rows), m_cols(cols) {
    allocate();
  }

  inline TinyMatrixXxX_(const TinyMatrixXxX_& other)
      : m_rows(other.m_rows), m_cols(other.m_cols) {
    allocate();
    for (int i = 0; i < m_cols; ++i) m_columns[i] = other.m_columns[i];
  }

  inline TinyMatrixXxX_& operator=(const TinyMatrixXxX_& other) {
    m_rows = other.m_rows;
    m_cols = other.m_cols;
    allocate();
    for (int i = 0; i < m_cols; ++i) m_columns[i] = other.m_columns[i];
    return *this;
  }


  virtual ~TinyMatrixXxX_() { delete[] m_columns; }

  void set_zero() {
    for (int i = 0; i < m_cols; ++i) m_columns[i].set_zero();
  }

  inline const TinyScalar& operator()(int row, int col) const {
    TinyConstants::FullAssert(0 <= row && row < m_rows);
    TinyConstants::FullAssert(0 <= col && col < m_cols);
    return m_columns[col][row];
  }
  inline TinyScalar& operator()(int row, int col) {
    TinyConstants::FullAssert(0 <= row && row < m_rows);
    TinyConstants::FullAssert(0 <= col && col < m_cols);
    return m_columns[col][row];
  }

  inline const ColumnVector& operator[](int col) const {
    TinyConstants::FullAssert(0 <= col && col < m_cols);
    return m_columns[col];
  }
  inline ColumnVector& operator[](int col) {
    TinyConstants::FullAssert(0 <= col && col < m_cols);
    return m_columns[col];
  }

  void print(const char* txt) const {
    printf("%s\n", txt);
    for (int r = 0; r < m_rows; r++) {
      for (int c = 0; c < m_cols; c++) {
        const TinyScalar& val = (*this)(r, c);

        TinyScalar v = TinyConstants::getDouble(val);
        printf("%f, ", v);
      }
      printf("\n");
    }
  }

  TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX> transpose() const {
    TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX> t(m_cols, m_rows);
    for (int i = 0; i < m_cols; ++i) {
      t.assign_vector_horizontal(i, 0, m_columns[i]);
    }
    return t;
  }


  TinyMatrixXxX_ operator-() {
    TinyMatrixXxX_ res(m_rows, m_cols);
    for (int j = 0; j < m_cols; ++j) {
      res.m_columns[j] = -m_columns[j];
    }
    return res;
  }

  TinyMatrixXxX_& operator+=(const TinyMatrixXxX_& m) {
    TinyConstants::FullAssert(m.m_cols == m_cols);
    TinyConstants::FullAssert(m.m_rows == m_rows);
    for (int j = 0; j < m_cols; ++j) {
      m_columns[j] += m.m_columns[j];
    }
    return *this;
  }

  TinyMatrixXxX_ operator+(const TinyMatrixXxX_& m) {
    TinyConstants::FullAssert(m.m_cols == m_cols);
    TinyConstants::FullAssert(m.m_rows == m_rows);
    TinyMatrixXxX_ res(m_rows, m_cols);
    for (int j = 0; j < m_cols; ++j) {
      res[j] = m_columns[j] + m.m_columns[j];
    }
    return res;
  }

  TinyMatrixXxX_ operator-(const TinyMatrixXxX_& m) {
    TinyConstants::FullAssert(m.m_cols == m_cols);
    TinyConstants::FullAssert(m.m_rows == m_rows);
    TinyMatrixXxX_ res(m_rows, m_cols);
    for (int j = 0; j < m_cols; ++j) {
      res[j] = m_columns[j] - m.m_columns[j];
    }
    return res;
  }
  /**
   * Multiples the LHS matrix times the RHS matrix.
   */
  friend TinyMatrixXxX_ operator*(const TinyMatrixXxX_& lhs,
                                  const TinyMatrixXxX_& rhs) {
    TinyConstants::FullAssert(lhs.m_cols == rhs.m_rows);
    TinyMatrixXxX_ res(lhs.m_rows, rhs.m_cols);
    res.set_zero();
    for (int i = 0; i < lhs.m_rows; ++i) {
      for (int j = 0; j < rhs.m_cols; ++j) {
        for (int k = 0; k < lhs.m_cols; ++k) {
          res(i, j) = res(i, j) + lhs(i, k) * rhs(k, j);
        }
      }
    }
    return res;
  }

  /**
   * Multiples the LHS matrix times the RHS vector.
   */
  template <template <typename, typename> typename VectorType>
  friend ColumnVector operator*(
      const TinyMatrixXxX_& lhs,
      const VectorType<TinyScalar, TinyConstants>& rhs) {
    TinyConstants::FullAssert(lhs.m_cols == rhs.m_size);
    ColumnVector res(lhs.m_rows);
    res.set_zero();

    for (int i = 0; i < lhs.m_rows; ++i) {
      for (int j = 0; j < lhs.m_cols; ++j) {
        res[i] = res[i] + lhs(i, j) * rhs[j];
      }
    }
    return res;
  }

  /**
   * Transposes the LHS matrix and multiplies it with the RHS vector.
   * RHS should be of matrix column dimension.
   */
  template <template <typename, typename> typename VectorType>
  TinyVectorX<TinyScalar, TinyConstants> mul_transpose(
      const VectorType<TinyScalar, TinyConstants>& rhs) {
    TinyConstants::FullAssert(m_rows == rhs.m_size);
    TinyVectorX<TinyScalar, TinyConstants> res(m_cols);
    for (int i = 0; i < m_cols; ++i) {
      res[i] = m_columns[i].dot(rhs);
    }
    return res;
  }

  // M(r,c) v(r) R(c)
  template <template <typename, typename> typename VectorType>
  void adj_mx_mul_transpose(
      const TinyVectorX<TinyScalar, TinyConstants>& R,
      const VectorType<TinyScalar, TinyConstants>& v,
      TinyMatrixXxX_<TinyScalar, TinyConstants, VectorType>& Rm,
      VectorType<TinyScalar, TinyConstants>& Rv) {
    Rm += vvt(v, R);
    Rv += (*this)*R;
  }


  // void adj_mx_mul_transpose(
  //     const TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVector3>& R,
  //     const TinyVector3<TinyScalar, TinyConstants>& v,
  //     TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVector3>& Rm,
  //     TinyVector3<TinyScalar, TinyConstants>& Rv) {
  //   Rm.m_columns += vvt(v, R);
  //   Rv += (*this)*R;
  // }


  // void adj_mx_mul_transpose(
  //     const TinyVectorX<TinyScalar, TinyConstants>& R,
  //     const TinyVectorX<TinyScalar, TinyConstants>& v,
  //     TinyMatrixXxX_& Rm,
  //     TinyVectorX<TinyScalar, TinyConstants>& Rv) {
  //   Rm.m_columns += vvt(v, R);
  //   Rv += (*this)*R;
  // }

  template <template <typename, typename> typename VectorType>
  void assign_vector_horizontal(
      int start_row_index, int start_col_index,
      const VectorType<TinyScalar, TinyConstants>& v) {
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index + v.m_size <= m_cols);
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index < m_rows);

    for (int i = 0; i < v.m_size; ++i) {
      (*this)(start_row_index, start_col_index + i) = v[i];
    }
  }

  template <template <typename, typename> typename VectorType>
  void assign_vector_horizontal_add(
      int start_row_index, int start_col_index,
      const VectorType<TinyScalar, TinyConstants>& v) {
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index + v.m_size <= m_cols);
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index < m_rows);

    for (int i = 0; i < v.m_size; ++i) {
      (*this)(start_row_index, start_col_index + i) += v[i];
    }
  }


  TinyVectorX<TinyScalar, TinyConstants> get_vector_horizontal(
      int start_row_index, int start_col_index, int length) const {
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index + length <= m_cols);
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index < m_rows);
    TinyVectorX<TinyScalar, TinyConstants> v(length);
    for (int i = 0; i < length; ++i) {
      v[i] = (*this)(start_row_index, i + start_col_index);
    }
    return v;
  }

  TinyVectorX<TinyScalar, TinyConstants> get_vector_vertical(
      int start_row_index, int start_col_index, int length) const {
    assert(start_row_index >= 0);
    assert(start_row_index + length <= m_rows);
    assert(start_col_index >= 0);
    assert(start_col_index < m_cols);
    TinyVectorX<TinyScalar, TinyConstants> v(length);
    for (int i = 0; i < length; ++i) {
      v[i] = (*this)(start_row_index + i,start_col_index);
    }
    return v;
  }

  template <template <typename, typename> typename VectorType>
  void assign_vector_vertical(int start_row_index, int start_col_index,
                              const VectorType<TinyScalar, TinyConstants>& v) {
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index + v.m_size <= m_rows);
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index < m_cols);

    ColumnVector& column = m_columns[start_col_index];
    for (int i = 0; i < v.m_size; ++i) {
      column[i + start_row_index] = v[i];
    }
  }

  template <typename MatrixType>
  void assign_matrix(int start_row_index, int start_col_index,
                     const MatrixType& m) {
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index + m.m_rows <= m_rows);
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index + m.m_cols <= m_cols);

    for (int j = 0; j < m.m_cols; ++j) {
      ColumnVector& column = m_columns[start_col_index + j];
      for (int i = 0; i < m.m_rows; ++i) {
        column[start_row_index + i] = m(i, j);
      }
    }
  }

  template <typename MatrixType>
  void assign_matrix_add(int start_row_index, int start_col_index,
                     const MatrixType& m) {
    TinyConstants::FullAssert(0 <= start_row_index);
    TinyConstants::FullAssert(start_row_index + m.m_rows <= m_rows);
    TinyConstants::FullAssert(0 <= start_col_index);
    TinyConstants::FullAssert(start_col_index + m.m_cols <= m_cols);

    for (int j = 0; j < m.m_cols; ++j) {
      ColumnVector& column = m_columns[start_col_index + j];
      for (int i = 0; i < m.m_rows; ++i) {
        column[start_row_index + i] += m(i, j);
      }
    }
  }

  TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX> block(
      int start_row_index, int start_col_index, int rows, int cols) const {
    assert(start_row_index >= 0);
    assert(start_row_index + rows <= m_rows);
    assert(start_col_index >= 0);
    assert(start_col_index + cols <= m_cols);
    TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX> m(rows, cols);
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j)
        m(i, j) = (*this)(i + start_row_index, j + start_col_index);
    }
    return m;
  }

  template <template <typename, typename> typename VectorType>
  static TinyMatrixXxX_<TinyScalar, TinyConstants, VectorType> vvt(
      const VectorType<TinyScalar, TinyConstants>& ma, 
      const TinyVectorX<TinyScalar, TinyConstants>& mb) {
    TinyMatrixXxX_<TinyScalar, TinyConstants, VectorType> 
              m(ma.m_size, mb.m_size);
    for (int i = 0; i < mb.m_size; i++) {
      m[i] = ma * mb[i];
    }

    return m;
  }

  /**
   *    main method for Cholesky decomposition.
   *    input/output  a  Symmetric positive def. matrix
   *    output        diagonal  vector of resulting diag of a
   *    inspired by public domain https://math.nist.gov/javanumerics/jama
   */
  bool cholesky_decomposition(
      TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& a,
      TinyVectorX<TinyScalar, TinyConstants>& diagonal) const {
    int i, j, k;
    TinyScalar sum;
    int n = a.m_cols;
    bool is_positive_definite = true;
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        sum = a[i][j];
        for (k = i - 1; k >= 0; k--) {
          sum = sum - a[i][k] * a[j][k];
        }
        if (i == j) {
          if (TinyConstants::getBool(
            sum <= TinyConstants::zero())) {
            is_positive_definite = false;
            break;
          }
          diagonal[i] = TinyConstants::sqrt1(sum);
        } else {
          a[j][i] = sum / diagonal[i];
        }
      }
    }
    return is_positive_definite;
  }

  /**
   *     Inverse of Cholesky decomposition.
   *
   *     input    A  Symmetric positive def. matrix
   *     output   a  inverse of lower decomposed matrix
   *     uses        cholesky_decomposition
   */
  bool inverse_cholesky_decomposition(
      const TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& A,
      TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& a) const {
    int i, j, k;
    int n = A.m_rows;
    TinyScalar sum;
    TinyVectorX<TinyScalar, TinyConstants> diagonal(A.m_rows);
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++) a[i][j] = A[i][j];
    bool is_positive_definite = cholesky_decomposition(a, diagonal);
    if (is_positive_definite) {
      for (i = 0; i < n; i++) {
        a[i][i] = TinyConstants::one() / diagonal[i];
        for (j = i + 1; j < n; j++) {
          sum = TinyConstants::zero();
          for (k = i; k < j; k++) {
            sum = sum - a[j][k] * a[k][i];
          }
          a[j][i] = sum / diagonal[j];
        }
      }
    }
    return is_positive_definite;
  }

  /**
   *     Inverse of a matrix, using Cholesky decomposition.
   *
   *     input    A  Symmetric positive def. matrix
   *     input    a  storage for the result
   *     output   boolean is_positive_definite if operation succeeded
   */
  bool inversed(
      TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& a) const {
    assert(m_cols == m_cols);
    assert(a.m_cols == m_cols);
    assert(a.m_rows == m_rows);

    const TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& A = *this;

    bool is_positive_definite = inverse_cholesky_decomposition(A, a);
    if (is_positive_definite) {
      int n = m_cols;
      int i, j, k;

      for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
          a[i][j] = TinyConstants::zero();
        }
      }

      for (i = 0; i < n; i++) {
        a[i][i] = a[i][i] * a[i][i];
        for (k = i + 1; k < n; k++) {
          a[i][i] = a[i][i] + a[k][i] * a[k][i];
        }
        for (j = i + 1; j < n; j++) {
          for (k = j; k < n; k++) {
            a[i][j] = a[i][j] + a[k][i] * a[k][j];
          }
        }
      }
      for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
          a[i][j] = a[j][i];
        }
      }
    }
    return is_positive_definite;
  }

  void adj_inversed(
      const TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& Rinv,
      const TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& b,
      TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>& R) const {

    R += - b.transpose() * Rinv * b.transpose();
  }
};

// template <typename TinyScalar, typename TinyConstants>
// using TinyMatrix6xX =
//     TinyMatrixXxX_<TinyScalar, TinyConstants, TinySpatialMotionVector>;

template <typename TinyScalar, typename TinyConstants>
using TinyMatrix3xX = TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVector3>;

template <typename TinyScalar, typename TinyConstants>
using TinyMatrixXxX = TinyMatrixXxX_<TinyScalar, TinyConstants, TinyVectorX>;

#endif  // TINY_MATRIXXXX_H
