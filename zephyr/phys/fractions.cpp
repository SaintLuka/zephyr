#include <iostream>
#include <zephyr/phys/fractions.h>

namespace zephyr::phys {

Fractions::Fractions(std::initializer_list<double> list) {
    if (list.size() > max_size) {
        throw std::runtime_error("When construct Fractions got list.size() > max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(max_size));
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

Fractions::Fractions(const std::array<double, max_size> &arr) : m_data(arr) {
    normalize();
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
    for (int i = 0; i < size(); ++i) {
        m_data[i] = 0.0;
    }
    if (0 <= idx && idx < size()) {
        m_data[idx] = 1.0;
    }
}

int Fractions::index() const {
    int idx = -1;
    for (int i = 0; i < size(); ++i) {
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

int Fractions::count() const {
    int idx = 0;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            idx = i;
        }
    }
    return idx + 1;
}

std::tuple<int, int> Fractions::pair() const {
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

std::ostream &operator<<(std::ostream &os, const Fractions &frac) {
    os << "{";
    for (int i = 0; i < frac.size() - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[frac.size()  - 1] << "}";
    return os;
}

ScalarSet::ScalarSet(std::initializer_list<double> list) {
    if (list.size() > max_size) {
        throw std::runtime_error("When construct ScalarSet got list.size() > max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(max_size));
    }

    int counter = 0;
    for (auto &elem: list) {
        m_data[counter++] = elem;
    }
    for (int i = counter; i < max_size; ++i) {
        m_data[i] = 0.0;
    }
}

ScalarSet::ScalarSet(const std::array<double, max_size>& arr)
    : m_data(arr) { }

ScalarSet::ScalarSet(const Fractions &frac)
    : m_data(frac.array()) { }

std::ostream &operator<<(std::ostream &os, const ScalarSet &frac) {
    os << "{";
    for (int i = 0; i < frac.size() - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[frac.size() - 1] << "}";
    return os;
}

VectorSet::VectorSet(std::initializer_list<Vector3d> list) {
    if (list.size() > max_size) {
        throw std::runtime_error("When construct VectorSet got list.size() > max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(max_size));
    }

    int counter = 0;
    for (auto &elem: list) {
        m_data[counter++] = elem;
    }
    for (int i = counter; i < max_size; ++i) {
        m_data[i] = Vector3d::Zero();
    }
}

std::ostream &operator<<(std::ostream &os, const VectorSet &arr) {
    os << "{ ";
    for (int i = 0; i < arr.size() - 1; ++i) {
        os << "{" << arr[i].transpose() << "}, ";
    }
    os << "{" << arr[arr.size() - 1].transpose() << "} }";
    return os;
}

} // namespace zephyr::phys