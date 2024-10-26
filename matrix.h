#ifndef MATRIX_H
#define MATRIX_H
#define MATRIX_SQUARE_MATRIX_IMPLEMENTED

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>

class MatrixOutOfRange : public std::exception {};
class MatrixIsDegenerateError : public std::exception {};

template <typename T, size_t X, size_t Y> class Matrix {
public:
  T matrix_[X][Y];

  [[nodiscard]] size_t RowsNumber() const { return X; }
  [[nodiscard]] size_t ColumnsNumber() const { return Y; }

  T &operator()(size_t x, size_t y) { return matrix_[x][y]; }
  const T &operator()(size_t x, size_t y) const { return matrix_[x][y]; }

  T &At(size_t x, size_t y) {
    if (0 <= x && x < X && 0 <= y && y < Y) {
      return this->matrix_[x][y];
    }
    throw MatrixOutOfRange{};
  }

  const T &At(size_t x, size_t y) const {
    if (0 <= x && x < X && 0 <= y && y < Y) {
      return this->matrix_[x][y];
    }
    throw MatrixOutOfRange{};
  }

  Matrix<T, X, Y> operator+(const Matrix<T, X, Y> &m) const {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = matrix_[i][j] + m(i, j);
      }
    }
    return new_m;
  }

  Matrix<T, X, Y> operator-(const Matrix<T, X, Y> &m) const {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = matrix_[i][j] - m(i, j);
      }
    }
    return new_m;
  }

  template <size_t Y2>
  Matrix<T, X, Y2> operator*(const Matrix<T, Y, Y2> &m) const {
    Matrix<T, X, Y2> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y2; j++) {
        T temp = T();
        for (size_t k = 0; k < Y; k++) {
          temp = temp + matrix_[i][k] * m(k, j);
        }
        new_m(i, j) = temp;
      }
    }
    return new_m;
  }

  Matrix<T, X, Y> &operator+=(const Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        matrix_[i][j] = matrix_[i][j] + m(i, j);
      }
    }
    return *this;
  }

  Matrix<T, X, Y> &operator-=(const Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        matrix_[i][j] = matrix_[i][j] - m(i, j);
      }
    }
    return *this;
  }

  Matrix<T, X, Y> &operator*=(const Matrix<T, Y, Y> &m) {
    Matrix<T, X, Y> new_m{};
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        T temp = T();
        for (size_t k = 0; k < Y; k++) {
          temp = temp + matrix_[i][k] * m(k, j);
        }
        new_m(i, j) = temp;
      }
    }

    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        (*this)(i, j) = new_m(i, j);
      }
    }
    return *this;
  }

  friend Matrix<T, X, Y> operator*(const Matrix<T, X, Y> &m, T num) {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = m(i, j) * num;
      }
    }
    return new_m;
  }

  friend Matrix<T, X, Y> operator/(const Matrix<T, X, Y> &m, T num) {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = m(i, j) / num;
      }
    }
    return new_m;
  }

  friend Matrix<T, X, Y> operator*(T num, const Matrix<T, X, Y> &m) {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = m(i, j) * num;
      }
    }
    return new_m;
  }

  friend Matrix<T, X, Y> operator/(T num, const Matrix<T, X, Y> &m) {
    Matrix<T, X, Y> new_m;
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        new_m(i, j) = m(i, j) / num;
      }
    }
    return new_m;
  }

  Matrix<T, X, Y> &operator*=(T num) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        matrix_[i][j] = matrix_[i][j] * num;
      }
    }
    return *this;
  }

  friend Matrix<T, X, Y> &operator*=(T num, const Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        m(i, j) = m(i, j) * num;
      }
    }
    return m;
  }

  friend Matrix<T, X, Y> &operator/=(T num, const Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        m(i, j) = m(i, j) / num;
      }
    }
    return m;
  }

  Matrix<T, X, Y> &operator/=(T num) {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        matrix_[i][j] = matrix_[i][j] / num;
      }
    }
    return *this;
  }

  bool operator==(const Matrix<T, X, Y> &m) const {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        if (matrix_[i][j] != m(i, j)) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const Matrix<T, X, Y> &m) const {
    for (size_t i = 0; i < X; i++) {
      for (size_t j = 0; j < Y; j++) {
        if (matrix_[i][j] != m(i, j)) {
          return true;
        }
      }
    }
    return false;
  }

  friend std::ostream &operator<<(std::ostream &os, const Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; ++i) {
      for (size_t j = 0; j < Y; ++j) {
        if (j != Y - 1) {
          os << m(i, j) << ' ';
        } else {
          os << m(i, j);
        }
      }
      os << std::endl;
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, Matrix<T, X, Y> &m) {
    for (size_t i = 0; i < X; ++i) {
      for (size_t j = 0; j < Y; ++j) {
        is >> m(i, j);
      }
    }
    return is;
  }
};

template <typename T, size_t X, size_t Y>
Matrix<T, Y, X> GetTransposed(const Matrix<T, X, Y> &m) {
  Matrix<T, Y, X> mt;
  for (size_t i = 0; i < X; i++) {
    for (size_t j = 0; j < Y; j++) {
      mt(j, i) = m(i, j);
    }
  }
  return mt;
}

template <typename T, size_t N>
void Transpose(Matrix<T, N, N> &matrix) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      std::swap(matrix(i, j), matrix(j, i));
    }
  }
}

template <typename T, size_t N>
T Trace(const Matrix<T, N, N> &matrix) {
  T s = 0;
  for (size_t i = 0; i < N; ++i) {
    s += matrix(i, i);
  }
  return s;
}

template <typename T, size_t N>
Matrix<T, N - 1, N - 1> GM(const Matrix<T, N, N> &matrix, size_t x, size_t y) {
  Matrix<T, N - 1, N - 1> m;
  for (size_t i = 0; i < N - 1; ++i) {
    for (size_t j = 0; j < N - 1; ++j) {
      size_t xj = (j >= x) ? j + 1 : j;
      size_t yi = (i >= y) ? i + 1 : i;
      m(i, j) = matrix(yi, xj);
    }
  }
  return m;
}


template <typename T, size_t N>
T Determinant(const Matrix<T, N, N> &matrix) {
  if (N == 0) {
    return 1;
  }
  Matrix<T, N, N> result = matrix;
  for (size_t k = 0; k < N - 1; ++k) {
    for (size_t i = k + 1; i < N; ++i) {
      for (size_t j = k + 1; j < N; ++j) {
        result(i, j) = result(k, k) * result(i, j) - result(i, k) * result(k, j);
        if (k != 0) {
          result(i, j) /= result(k - 1, k - 1);
        }
      }
    }
  }
  return result(N - 1, N - 1);
}

template <typename T, size_t N>
Matrix<T, N, N> GetInversed(const Matrix<T, N, N> &matrix) {
  T det = Determinant(matrix);
  if (det == 0) {
    throw MatrixIsDegenerateError{};
  }
  Matrix<T, N, N> answ;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      T zn = ((i + j) % 2 == 0) ? 1 : -1;

      answ(i, j) = zn * Determinant(GM(matrix, j, i));
    }
  }
  Transpose(answ);
  answ /= det;
  return answ;
}
template <typename T>
Matrix<T, 1, 1> GetInversed(const Matrix<T, 1, 1> &matrix) {
  T det = Determinant(matrix);
  if (det == 0) {
    throw MatrixIsDegenerateError{};
  }
  return matrix / det / det;
}

template <typename T, size_t N>
void Inverse(Matrix<T, N, N> &matrix) {
  T det = Determinant(matrix);
  if (det == 0) {
    throw MatrixIsDegenerateError{};
  }
  Matrix<T, N, N> answ;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      T zn = ((i + j) % 2 == 0) ? 1 : -1;

      
      answ(i, j) = zn * Determinant(GM(matrix, j, i));
    }
  }
  Transpose(answ);
  answ /= det;
  matrix = answ;
}

template <typename T>
void Inverse(Matrix<T, 1, 1> &matrix) {
  T det = Determinant(matrix);
  if (det == 0) {
    throw MatrixIsDegenerateError{};
  }
  matrix /= det;
  matrix /= det;
}

#endif // MATRIX_H
