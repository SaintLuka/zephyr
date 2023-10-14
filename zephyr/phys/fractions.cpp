#include <iostream>
#include <zephyr/phys/fractions.h>

namespace zephyr::phys {

Fractions::Fractions() {
    m_data.fill(0.0);
}

Fractions::Fractions(std::initializer_list<double> list) : m_begin_size(list.size()) {
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

Fractions::Fractions(const std::vector<double> &list) : m_begin_size(list.size()) {
    if (list.size() > Fractions::max_size) {
        throw std::runtime_error("When construct Fractions got list.size() > Fractions::max_size (" +
                                 std::to_string(list.size()) + " > " + std::to_string(Fractions::max_size));
    }

    for (size_t i = 0; i < list.size(); ++i) {
        m_data[i] = list[i];
    }
    for (int i = list.size(); i < max_size; ++i) {
        m_data[i] = 0.0;
    }

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
    for (int i = 0; i < m_begin_size; ++i) {
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
    for (int i = 0; i < m_begin_size; ++i) {
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
    double sum = 0.0;
    for (int i = 0; i < m_begin_size; ++i) {
        sum += m_data[i];
    }
    for (int i = 0; i < m_begin_size; ++i) {
        m_data[i] /= sum;
    }
}

void Fractions::cutoff(double eps) {
    double sum = 0.0;
    for (int i = 0; i < m_begin_size; ++i) {
        if (m_data[i] < eps) {
            m_data[i] = 0.0;
        } else if (m_data[i] > 1.0 - eps) {
            m_data[i] = 1.0;
        }
        sum += m_data[i];
    }
    for (int i = 0; i < m_begin_size; ++i) {
        m_data[i] /= sum;
    }
}

size_t Fractions::get_begin_size() const {
    return m_begin_size;
}

std::ostream &operator<<(std::ostream &os, const Fractions &frac) {
    os << "[ ";
    for (int i = 0; i < frac.get_begin_size() - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[frac.get_begin_size() - 1] << " ]";
    return os;
}

} // namespace zephyr