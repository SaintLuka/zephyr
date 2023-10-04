#include <iostream>
#include <zephyr/phys/fractions.h>

namespace zephyr { namespace phys {

Fractions::Fractions() {
    m_data.fill(0.0);
}

Fractions::Fractions(std::initializer_list<double> list) {
    int counter = 0;
    for (auto elem: list) {
        if (counter < max_size) {
            m_data[counter] = elem;
        }
        else {
            break;
        }
        ++counter;
    }
    for (int i = counter; i < max_size; ++i) {
        m_data[i] = 0.0;
    }
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

bool Fractions::is_pure() const {
    bool pure = false;
    for (int i = 0; i < max_size; ++i) {
        if (has(i)) {
            if (!pure) {
                // Нашли первое вещество
                pure = true;
            }
            else {
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
            }
            else {
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
    for (int i = 0; i < max_size; ++i) {
        sum += m_data[i];
    }
    for (int i = 0; i < max_size; ++i) {
        m_data[i] /= sum;
    }
}

void Fractions::cutoff(double eps) {
    double sum = 0.0;
    for (int i = 0; i < max_size; ++i) {
        if (m_data[i] < eps) {
            m_data[i] = 0.0;
        }
        else if (m_data[i] > 1.0 - eps) {
            m_data[i] = 1.0;
        }
        sum += m_data[i];
    }
    for (int i = 0; i < max_size; ++i) {
        m_data[i] /= sum;
    }
}

std::ostream& operator<<(std::ostream& os, const Fractions& frac) {
    os << "[ ";
    for (int i = 0; i < Fractions::max_size - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[Fractions::max_size - 1] << " ]";
    return os;
}

} // namespace phys
} // namespace zephyr