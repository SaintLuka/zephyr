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

Fractions::Fractions(const FractionsFlux &frac_flux) {
    std::copy(frac_flux.m_data.begin(), frac_flux.m_data.end(), m_data.begin());

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
    double sub = 0.0;
    for (auto &v: m_data) {
        sub += abs(v);
    }

    double sum = 0;
    for (auto &v: m_data) {
        if (v < 0) {
            if (-1e-5 * sub < v)
                v = 0;
            else {
                std::cerr << "mass_frac[i] < 0: " << *this;
//                exit(1);
                throw std::runtime_error("mass_frac[i] < 0");
            }
        } else
            sum += v;
    }

    for (int i = 0; i < max_size; ++i) {
        m_data[i] /= sum;
    }
    cutoff(1e-12);
}

void Fractions::cutoff(double eps) {
    double sum = 0.0;
    for (int i = 0; i < max_size; ++i) {
        if (m_data[i] < eps) {
            m_data[i] = 0.0;
        } else if (m_data[i] > 1.0 - eps) {
            m_data[i] = 1.0;
        }
        sum += m_data[i];
    }
    for (int i = 0; i < max_size; ++i) {
        m_data[i] /= sum;
    }
}

bool Fractions::empty() const {
    for (auto &v: m_data)
        if (v > 1e-8)
            return false;

    return true;
}

size_t Fractions::get_size() const {
    return max_size;
}

void Fractions::fix() {
    for (auto &v: m_data)
        if (v < 0)
            v = 0;
        else if (v > 1)
            v = 1;

    normalize();
}

std::ostream &operator<<(std::ostream &os, const Fractions &frac) {
    os << "{";
    for (int i = 0; i < frac.get_size() - 1; ++i) {
        os << frac[i] << ", ";
    }
    os << frac[frac.get_size() - 1] << "}";
    return os;
}

FractionsFlux::FractionsFlux() {
    m_data.fill(0);
}

FractionsFlux::FractionsFlux(const Fractions &frac) : m_data(frac.get_data()) {}

FractionsFlux::FractionsFlux(const std::vector<double> &vec) {
    if (vec.size() > Fractions::max_size) {
        throw std::runtime_error("When construct FractionsFlux got vec.size() > Fractions::max_size (" +
                                 std::to_string(vec.size()) + " > " + std::to_string(Fractions::max_size));
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        m_data[i] = vec[i];
    }
    for (int i = vec.size(); i < Fractions::max_size; ++i) {
        m_data[i] = 0.0;
    }
}

bool FractionsFlux::has(int idx) const {
    if (idx > m_data.size())
        return false;
    return m_data[idx] > 0.0;
}

double &FractionsFlux::operator[](size_t idx) {
    return m_data[idx];
}

const double &FractionsFlux::operator[](size_t idx) const {
    return m_data[idx];
}

std::ostream &operator<<(std::ostream &os, const FractionsFlux &frac) {
    os << "{";
    for (int i = 0; i < frac.get_size() - 1; ++i) {
        os << frac.m_data[i] << ", ";
    }
    os << frac.m_data[frac.get_size() - 1] << "}";
    return os;
}

} // namespace zephyr