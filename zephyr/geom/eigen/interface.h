/// @brief Плохая копия интерфейса Eigen, потому что Eigen иногда тупит,
/// выдает странные ошибки, которые тяжело отлаживаются.
#pragma once

#include <iostream>
#include <array>
#include <cassert>

#include <Eigen/Dense>

namespace zephyr::geom::eigen_wrapper {

template <typename _Scalar, int _Size, int _Unused = _Size>
class DiagonalMatrix;

/// @brief Класс матрицы. Вектор при _Cols = 1
template<typename _Scalar, int _Rows, int _Cols, int _Size = _Rows * _Cols,
        bool _IsVector = (_Rows == 1 || _Cols == 1),
        bool _IsVector3d = ((_Rows == 3 || _Cols == 3) && _Size == 3)>
class Matrix {
public:
    using SameMatrix = Matrix<_Scalar, _Rows, _Cols, _Size>;

    /// @brief Инициализирую нулями, хотя Eigen по умолчанию
    /// так не делает!!
    Matrix() {
        for (int i = 0; i < _Size; ++i) {
            get(i) = 0.0;
        }
    }

    /// @brief Конструктор векторов до четырехмерных
    Matrix(const _Scalar& _x,
           const _Scalar& _y = 0.0,
           const _Scalar& _z = 0.0,
           const _Scalar& _w = 0.0) : Matrix() {
        x() = _x;
        if (_Size > 1) {
            y() = _y;
            if (_Size > 2) {
                z() = _z;
                if (_Size > 3) {
                    w() = _w;
                }
            }
        }
    }

    /// @brief Возвращает нулевую матрицу
    static SameMatrix Zero() {
        SameMatrix sm;
        for (int i = 0; i < _Size; ++i) {
            sm.get(i) = 0.0;
        }
        return sm;
    }

    /// @brief Возвращает нулевую матрицу
    static SameMatrix Identity() {
        SameMatrix sm;
        for (int i = 0; i < std::min(_Rows, _Cols); ++i) {
            sm(i, i) = 1.0;
        }
        return sm;
    }

    static SameMatrix UnitX() {
        return SameMatrix(1.0);
    }

    static SameMatrix UnitY() {
        return SameMatrix(0.0, 1.0);
    }

    static SameMatrix UnitZ() {
        return SameMatrix(0.0, 0.0, 1.0);
    }

    /// @brief Указатель на данные
    inline _Scalar* data() {
        return m_data[0].data();
    }

    /// @brief Константный указатель на данные
    inline const _Scalar* data() const {
        return m_data[0].data();
    }

    /// @brief Данные по одномерному индексу
    inline _Scalar& get(int idx) {
        assert(idx < _Size);
        return data()[idx];
    }

    /// @brief Данные по одномерному индексу
    inline const _Scalar& get(int idx) const {
        assert(idx < _Size);
        return data()[idx];
    }

    /// @brief Оператор доступа по единственному индексу [].
    /// Актуально для векторов.
    inline _Scalar& operator[](int i) {
        if (!_IsVector) {
            throw std::runtime_error("eigen interface bad []");
        }
        return get(i);
    }

    /// @brief Оператор доступа по единственному индексу [].
    /// Актуально для векторов.
    inline const _Scalar& operator[](int i) const {
        if (!_IsVector) {
            throw std::runtime_error("eigen interface bad []");
        }
        return get(i);
    }

    /// @brief Оператор доступа по паре индексов (i, j)
    inline _Scalar& operator()(int row, int col) {
        assert(row < _Rows);
        assert(col < _Cols);
        return m_data[row][col];
    }

    /// @brief Оператор доступа по паре индексов (i, j)
    inline const _Scalar& operator()(int row, int col) const {
        assert(row < _Rows);
        assert(col < _Cols);
        return m_data[row][col];
    }

    inline _Scalar& x() { return get(0); }
    inline _Scalar& y() { return get(1); }
    inline _Scalar& z() { return get(2); }
    inline _Scalar& w() { return get(3); }

    inline const _Scalar& x() const { return get(0); }
    inline const _Scalar& y() const { return get(1); }
    inline const _Scalar& z() const { return get(2); }
    inline const _Scalar& w() const { return get(3); }


    SameMatrix operator-() const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = -get(i);
        }
        return res;
    }

    SameMatrix operator+(const SameMatrix& r) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = get(i) + r.get(i);
        }
        return res;
    }

    SameMatrix operator-(const SameMatrix& r) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = get(i) - r.get(i);
        }
        return res;
    }

    SameMatrix operator*(const _Scalar& q) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = q * get(i);
        }
        return res;
    }

    SameMatrix operator/(const _Scalar& q) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = get(i) / q;
        }
        return res;
    }

    const SameMatrix& operator+=(const SameMatrix& r) {
        for (int i = 0; i < _Size; ++i) {
            get(i) += r.get(i);
        }
        return *this;
    }

    const SameMatrix& operator-=(const SameMatrix& r) {
        for (int i = 0; i < _Size; ++i) {
            get(i) -= r.get(i);
        }
        return *this;
    }

    const SameMatrix& operator*=(const _Scalar& q) {
        for (int i = 0; i < _Size; ++i) {
            get(i) *= q;
        }
        return *this;
    }

    const SameMatrix& operator/=(const _Scalar& q) {
        for (int i = 0; i < _Size; ++i) {
            get(i) /= q;
        }
        return *this;
    }

    template <int N>
    Matrix<_Scalar, _Rows, N> operator*(const Matrix<_Scalar, _Cols, N>& r) const {
        Matrix<_Scalar, _Rows, N> res;
        for (int i = 0; i < _Rows; ++i) {
            for (int j = 0; j < N; ++j) {
                res(i, j) = 0.0;
                for (int k = 0; k < _Cols; ++k) {
                    res(i, j) += m_data[i][k] * r(k, j);
                }
            }
        }
        return res;
    }

    SameMatrix operator*(const DiagonalMatrix<_Scalar, _Cols>& diag) const {
        SameMatrix res;
        for (int i = 0; i < _Rows; ++i) {
            for (int j = 0; j < _Cols; ++j) {
                res(i, j) = m_data[i][j] * diag[j];
            }
        }
        return res;
    }

    _Scalar squaredNorm() const {
        _Scalar res(0);
        for (int i = 0; i < _Size; ++i) {
            res += get(i) * get(i);
        }
        return res;
    }

    _Scalar norm() const {
        return std::sqrt(squaredNorm());
    }

    void normalize() {
        _Scalar sn = squaredNorm();
        if (sn > _Scalar(0)) {
            *this /= std::sqrt(sn);
        }
    }

    SameMatrix normalized() const {
        _Scalar sn = squaredNorm();
        if (sn > _Scalar(0)) {
            return *this / std::sqrt(sn);
        }
        return *this;
    }

    _Scalar dot(const SameMatrix& r) const {
        static_assert(_IsVector, "eigen interface bad dot");

        _Scalar res(0);
        for (int i = 0; i < _Size; ++i) {
            res += get(i) * r.get(i);
        }
        return res;
    }

    SameMatrix cross(const SameMatrix& b) const {
        static_assert(_IsVector3d, "eigen interface bad cross");

        const auto& a = *this;
        SameMatrix c = {
                a.y() * b.z() - b.y() * a.z(),
                b.x() * a.z() - a.x() * b.z(),
                a.x() * b.y() - b.x() * a.y()
        };
        return c;
    }

    SameMatrix inverse() const {
        Eigen::Matrix<double, _Rows, _Cols> mat;
        for (int i = 0; i < _Size; ++i) {
            mat.data()[i] = get(i);
        }
        mat = mat.inverse();

        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = mat.data()[i];
        }
        return res;
    }

    bool hasNaN() const {
        for (int i = 0; i < _Size; ++i) {
            if (std::isnan(get(i))) {
                return true;
            }
        }
        return false;
    }

    bool isZero() const {
        for (int i = 0; i < _Size; ++i) {
            if (get(i) != 0.0) {
                return false;
            }
        }
        return true;
    }

    SameMatrix cwiseAbs() const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = std::abs(get(i));
        }
        return res;
    }

    SameMatrix cwiseMin(const SameMatrix& m) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = std::min(get(i), m.get(i));
        }
        return res;
    }

    SameMatrix cwiseMax(const SameMatrix& m) const {
        SameMatrix res;
        for (int i = 0; i < _Size; ++i) {
            res.get(i) = std::max(get(i), m.get(i));
        }
        return res;
    }

    DiagonalMatrix<_Scalar, _Size> asDiagonal() const {
        DiagonalMatrix<_Scalar, _Size> diag;
        for (int i = 0; i < _Size; ++i) {
            diag[i] = get(i);
        }
        return diag;
    }

    Matrix<_Scalar, _Cols, _Rows> transpose() const {
        Matrix<_Scalar, _Cols, _Rows> res;
        for (int i = 0; i < _Rows; ++i) {
            for (int j = 0; j < _Cols; ++j) {
                res(j, i) = m_data[i][j];
            }
        }
        return res;
    }

    struct CommaInput {
    public:
        CommaInput(_Scalar* _ptr, int _rest)
                : ptr(_ptr), rest(_rest) {

        }

        CommaInput& operator,(_Scalar x) {
            if (rest > 0) {
                *(ptr++) = x;
                --rest;
            }
            return *this;
        }

        _Scalar* ptr;
        int rest;
    };

    CommaInput operator<<(double x) {
        CommaInput ci(data(), _Size);
        return ci, x;
    }

private:

    /// @brief Двумерный массив с данными, для удобства отладки
    /// в Eigen данные хранятся в однмерном массиве
    std::array<std::array<_Scalar, _Cols>, _Rows> m_data;
};

template<typename _Scalar, int _Rows, int _Cols>
std::ostream& operator<<(std::ostream& os, const Matrix<_Scalar, _Rows, _Cols>& m) {
    int _Size = _Rows * _Cols;
    if (_Rows == 1 || _Cols == 1) {
        for (int i = 0; i < _Size - 1; ++i) {
            os << m[i];
            if (i < _Size - 1) {
                os << " ";
            }
        }
        os << m[_Size - 1];
    }
    else {
        for (int i = 0; i < _Rows; ++i) {
            for (int j = 0; j < _Cols; ++j) {
                os << m(i, j);
                if (j < _Cols - 1) {
                    os << " ";
                }
            }
            if (i < _Rows - 1) {
                os << "\n";
            }
        }
    }
    return os;
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar, _Rows, _Cols> operator+(
        const Matrix<_Scalar, _Rows, _Cols>& m1,
        const Matrix<_Scalar, _Rows, _Cols>& m2) {
    Matrix<_Scalar, _Rows, _Cols> res;
    for (int i = 0; i < _Rows; ++i) {
        for (int j = 0; j < _Cols; ++j) {
            res(i, j) = m1(i, j) + m2(i, j);
        }
    }
    return res;
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar, _Rows, _Cols> operator-(
        const Matrix<_Scalar, _Rows, _Cols>& m1,
        const Matrix<_Scalar, _Rows, _Cols>& m2) {
    Matrix<_Scalar, _Rows, _Cols> res;
    for (int i = 0; i < _Rows; ++i) {
        for (int j = 0; j < _Cols; ++j) {
            res(i, j) = m1(i, j) - m2(i, j);
        }
    }
    return res;
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar, _Rows, _Cols> operator*(
        const _Scalar& q, const Matrix<_Scalar, _Rows, _Cols>& m) {
    Matrix<_Scalar, _Rows, _Cols> res;
    for (int i = 0; i < _Rows; ++i) {
        for (int j = 0; j < _Cols; ++j) {
            res(i, j) = q * m(i, j);
        }
    }
    return res;
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar, _Rows, _Cols> operator*(
        int q, const Matrix<_Scalar, _Rows, _Cols>& m) {
    return _Scalar(q) * m;
}

template <typename _Scalar, int _Size, int _Unused>
class DiagonalMatrix {
public:
    /// @brief Инициализирую нулями, хотя Eigen по умолчанию
    /// так не делает!!
    DiagonalMatrix() {
        for (int i = 0; i < _Size; ++i) {
            get(i) = 0.0;
        }
    }

    /// @brief Возвращает нулевую матрицу
    static DiagonalMatrix Zero() {
        DiagonalMatrix sm;
        for (int i = 0; i < _Size; ++i) {
            sm.get(i) = 0.0;
        }
        return sm;
    }

    /// @brief Указатель на данные
    inline _Scalar* data() {
        return m_data.data();
    }

    /// @brief Константный указатель на данные
    inline const _Scalar* data() const {
        return m_data.data();
    }

    /// @brief Оператор доступа по единственному индексу [].
    /// Актуально для векторов.
    inline _Scalar& operator[](int i) {
        return get(i);
    }

    /// @brief Оператор доступа по единственному индексу [].
    /// Актуально для векторов.
    inline const _Scalar& operator[](int i) const {
        return get(i);
    }

    /// @brief Оператор доступа по единственному индексу ().
    /// Актуально для векторов.
    inline _Scalar& operator()(int i) {
        return get(i);
    }

    /// @brief Оператор доступа по единственному индексу ().
    /// Актуально для векторов.
    inline const _Scalar& operator()(int i) const {
        return get(i);
    }

    template <int N>
    Matrix<_Scalar, _Size, N> operator*(const Matrix<_Scalar, _Size, N>& r) const {
        Matrix<_Scalar, _Size, N> res;
        for (int i = 0; i < _Size; ++i) {
            for (int j = 0; j < N; ++j) {
                res(i, j) = 0.0;
                for (int k = 0; k < _Size; ++k) {
                    res(i, j) += get(i) * r(k, j);
                }
            }
        }
        return res;
    }


private:
    /// @brief Данные по одномерному индексу
    inline _Scalar& get(int idx) {
        assert(idx < _Size);
        return data()[idx];
    }

    /// @brief Данные по одномерному индексу
    inline const _Scalar& get(int idx) const {
        assert(idx < _Size);
        return data()[idx];
    }

    /// @brief Данные
    std::array<_Scalar, _Size> m_data;
};

template <typename _Scalar, int _Size, int _Unused>
using Array = Matrix<_Scalar, _Size, _Unused>;

} // namespace zephyr::geom