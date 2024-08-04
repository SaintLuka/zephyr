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

void Fractions::normalize() {
    cutoff(1e-14);
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
        if (v > 1e-14) {
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

double &ScalarSet::operator[](size_t idx) {
    return m_data[idx];
}

const double &ScalarSet::operator[](size_t idx) const {
    return m_data[idx];
}

std::ostream &operator<<(std::ostream &os, const ScalarSet &frac) {
    os << "{";
    for (int i = 0; i < Fractions::size() - 1; ++i) {
        os << frac.m_data[i] << ", ";
    }
    os << frac.m_data[Fractions::size() - 1] << "}";
    return os;
}

} // namespace zephyr::phys