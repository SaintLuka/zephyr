#include <iostream>
#include <zephyr/phys/fractions.h>

namespace zephyr::phys {

Fractions::Fractions() {
    m_data.fill(0.0);
}

Fractions::Fractions(std::initializer_list<double> list) {
    if (list.size() > Fractions::max_size) {
        throw std::runtime_error("When construct Fractions got list.size() > Fractions::max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(Fractions::max_size));
    }

    int counter = 0;
    for (auto &elem: list) {
        m_data[counter++] = elem;
    }
    for (int i = counter; i < max_size; ++i) {
        m_data[i] = 0.0;
    }

    normalize();
}

Fractions::Fractions(const std::vector<double> &vec) {
    if (vec.size() > Fractions::max_size) {
        throw std::runtime_error("When construct Fractions got vec.size() > Fractions::max_size (" +
                                 std::to_string(vec.size()) + " > " + std::to_string(Fractions::max_size));
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        m_data[i] = vec[i];
    }
    for (int i = vec.size(); i < max_size; ++i) {
        m_data[i] = 0.0;
    }

    normalize();
}

Fractions::Fractions(const std::array<double, max_size> &arr) : m_data(arr) {
    normalize();
}

Fractions::Fractions(const ScalarSet &scalars) {
    std::copy(scalars.m_data.begin(), scalars.m_data.end(), m_data.begin());

    normalize();
}

bool Fractions::has(int idx) const {
    return m_data[idx] > 0.0;
}

double &Fractions::operator[](int idx) {
    return m_data[idx];
}

const double &Fractions::operator[](int idx) const {
    return m_data[idx];
}

double &Fractions::operator[](size_t idx) {
    return m_data[idx];
}

const double &Fractions::operator[](size_t idx) const {
    return m_data[idx];
}

bool Fractions::is_pure() const {
    bool pure = false;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            if (!pure) {
                // Нашли первое вещество
                pure = true;
            } else {
                // Нашли второе вещество
                return false;
            }
        }
    }
    // Либо не нашли вещество (тогда false),
    // Либо нашли только одно, тогда true
    return pure;
}

void Fractions::set_pure(int idx) {
    for (int i = 0; i < Fractions::max_size; ++i) {
        m_data[i] = 0.0;
    }
    if (0 <= idx && idx < Fractions::max_size) {
        m_data[idx] = 1.0;
    }
}

int Fractions::index() const {
    int idx = -1;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            if (idx < 0) {
                // Нашли первое вещество
                idx = i;
            } else {
                // Нашли второе вещество
                return -1;
            }
        }
    }
    // Либо не нашли вещество (тогда -1),
    // Либо нашли только одно, тогда i
    return idx;
}

int Fractions::count() const {
    int idx = 0;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            idx = i;
        }
    }
    return idx + 1;
}

int Fractions::nonzero() const {
    int counter = 0;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            ++counter;
        }
    }
    return counter;
}

std::array<int, 2> Fractions::pair() const {
    int first = -1;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            if (first < 0) {
                // Нашли первое вещество
                first = i;
            } else {
                // Нашли второе вещество
                return {first, i};
            }
        }
    }
    // Либо не нашли вещество,
    // либо нашли только одно
    return {-1, -1};
}

void Fractions::normalize() {
    cutoff(minimal);
}

void Fractions::cutoff(double eps) {
    double sum = 0.0;
    for (int i = 0; i < max_size; ++i) {
        if (m_data[i] < eps) {
            m_data[i] = 0.0;
        } else if (m_data[i] > 1.0) {
            m_data[i] = 1.0;
        }
        sum += m_data[i];
    }
    for (int i = 0; i < max_size; ++i) {
        m_data[i] /= sum;
    }
}

bool Fractions::empty() const {
    for (auto &v: m_data) {
        if (v > minimal) {
            return false;
        }
    }

    return true;
}

std::ostream &operator<<(std::ostream &os, const Fractions &frac) {
    os << "{";
    for (int i = 0; i < Fractions::size() - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[Fractions::size()  - 1] << "}";
    return os;
}

ScalarSet::ScalarSet() {
    m_data.fill(0.0);
}

ScalarSet::ScalarSet(double val) {
    m_data.fill(val);
}

ScalarSet::ScalarSet(std::initializer_list<double> list) {
    if (list.size() > Fractions::max_size) {
        throw std::runtime_error("When construct ScalarSet got list.size() > Fractions::max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(Fractions::max_size));
    }

    int counter = 0;
    for (auto &elem: list) {
        m_data[counter++] = elem;
    }
    for (int i = counter; i < Fractions::max_size; ++i) {
        m_data[i] = 0.0;
    }
}

ScalarSet::ScalarSet(double val, int idx) {
    m_data.fill(0.0);
    m_data[idx] = val;
}

ScalarSet::ScalarSet(const Fractions &frac)
    : m_data(frac.data_ref()) {}

ScalarSet::ScalarSet(const std::vector<double> &vec) {
    if (vec.size() > Fractions::max_size) {
        throw std::runtime_error("When construct ScalarSet got vec.size() > Fractions::max_size (" +
                                 std::to_string(vec.size()) + " > " + std::to_string(Fractions::max_size));
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        m_data[i] = vec[i];
    }
    for (int i = vec.size(); i < Fractions::max_size; ++i) {
        m_data[i] = 0.0;
    }
}

std::ostream &operator<<(std::ostream &os, const ScalarSet &frac) {
    os << "{";
    for (int i = 0; i < Fractions::size() - 1; ++i) {
        os << frac.m_data[i] << ", ";
    }
    os << frac.m_data[Fractions::size() - 1] << "}";
    return os;
}

VectorSet::VectorSet() {
    m_data.fill(Vector3d::Zero());
}

VectorSet::VectorSet(const Vector3d& vec, int idx) {
    m_data.fill(Vector3d::Zero());
    m_data[idx] = vec;
}

VectorSet::VectorSet(const std::vector<Vector3d> &list) {
    if (list.size() > Fractions::max_size) {
        throw std::runtime_error("When construct VectorSet got list.size() > Fractions::max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(Fractions::max_size));
    }

    int counter = 0;
    for (auto &elem: list) {
        m_data[counter++] = elem;
    }
    for (int i = counter; i < Fractions::max_size; ++i) {
        m_data[i] = Vector3d::Zero();
    }
}

std::ostream &operator<<(std::ostream &os, const VectorSet &frac) {
    os << "{ ";
    for (int i = 0; i < Fractions::size() - 1; ++i) {
        os << "{" << frac.m_data[i].transpose() << "}, ";
    }
    os << "{" << frac.m_data[Fractions::size() - 1].transpose() << "} }";
    return os;
}

} // namespace zephyr::phys