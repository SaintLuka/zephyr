#include <set>
#include <iostream>
#include <iomanip>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <zephyr/math/cfd/limiter.h>
#include <cmath>

namespace zephyr { namespace math {

namespace limiters {

inline double sqr(double x) {
    return x * x;
}

double zero(double r) {
    return 0.0;
}

double zero(double a, double b) {
    return 0.0;
}

double one(double r) {
    return 1.0;
}

double one(double a, double b) {
    return 1.0;
}

double minmod(double r) {
    if (std::isnan(r)) { return 0.0; }
    return std::max(0.0, std::min(r, 1.0));
}

double minmod(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    return std::min(a / b, 1.0);
}

double MC(double r) {
    if (std::isnan(r)) { return 0.0; }
    return std::max(0.0, std::min(2 * r, std::min(0.5 * (1.0 + r), 2.0)));
}

double MC(double a, double b) {
    return MC(a / b);
}

double vanLeer(double r) {
    if (r <= 0.0) { return 0.0; }
    if (std::isnan(r)) { return 0.0; }
    if (std::isinf(r)) { return 2.0; }

    return 2.0 * r / (r + 1.0);
}

double vanLeer(double a, double b) {
    if (a * b < 0.0) { return 0.0; }
    if (a == 0.0) { return 0.0; }

    return 2 * std::abs(a) / (std::abs(a) + std::abs(b));
}

double vanAlbada1(double r) {
    if (r <= 0.0) { return 0.0; }
    if (std::isnan(r)) { return 0.0; }
    if (std::isinf(r)) { return 1.0; }

    return r * (r + 1.0) / (sqr(r) + 1.0);
}

double vanAlbada1(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    if (b == 0.0) { return 1.0; }

    return a * (a + b) / (sqr(a) + sqr(b));
}

double vanAlbada2(double r) {
    if (r <= 0.0) { return 0.0; }
    if (std::isnan(r)) { return 0.0; }
    if (std::isinf(r)) { return 0.0; }

    return 2.0 * r / (sqr(r) + 1.0);
}

double vanAlbada2(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    double eps = 10e-6;
    return (2.0 * a * b + eps) / (sqr(a) + sqr(b) + eps);
}

double CHARM(double r) {
    if (std::isnan(r)) { return 0.0; }
    return r > 0.0 ? (3.0 + 1.0 / r) / sqr(r + 1.0 / r) : 0.0;
}

double CHARM(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    return a * (3 * a + b) / sqr(a + b);
}

double HCUS(double r) {
    if (std::isnan(r)) { return 0.0; }
    return r > 0.0 ? 3.0 / (1.0 + 2.0 / r) : 0.0;
}

double HCUS(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    return (3.0 * a) / (a + 2.0 * b);
}

double HQUICK(double r) {
    if (std::isnan(r)) { return 0.0; }
    return r > 0.0 ? 4.0 / (1.0 + 3.0 / r) : 0.0;
}

double HQUICK(double a, double b) {
    if (a * b <= 0.0) { return 0.0; }
    return (4.0 * a) / (a + 3.0 * b);
}

double Koren(double r) {
    if (std::isnan(r)) { return 0.0; }
    return std::max(0.0, std::min(2 * r, std::min(2.0, (1 + 2 * r) / 3.0)));
}

double Koren(double a, double b) {
    return Koren(a / b);
}

double smart(double r) {
    if (std::isnan(r)) { return 0.0; }
    return std::max(0.0, std::min(2.0 * r, std::min(0.25 + 0.75 * r, 4.0)));
}

double smart(double a, double b) {
    return smart(a / b);
}

double superbee(double r) {
    if (std::isnan(r)) { return 0.0; }
    return std::max(0.0, std::min(std::min(2 * r, 1.0), std::min(r, 2.0)));
}

double superbee(double a, double b) {
    return superbee(a / b);
}

} // namespace limiters

Limiter::Limiter()
    : Limiter("minmod") { }

Limiter::Limiter(const char *name)
    : Limiter(std::string(name)) { }

Limiter::Limiter(const std::string& orig_name) {
    using namespace boost::algorithm;

    std::string name = orig_name;
    to_lower(name);
    erase_all(name, " ");
    erase_all(name, "-");
    erase_all(name, "_");
    erase_all(name, "lim");
    erase_all(name, "limiter");

    m_active = true;

    if (name.empty() ||
        contains(name, "off") ||
        contains(name, "none") ||
        contains(name, "zero") ||
        contains(name, "empty")) {
        m_name = "none";
        m_active = false;
        func1 = limiters::zero;
        func2 = limiters::zero;
    }
    else if (contains(name, "one") ||
             contains(name, "unit")) {
        m_name = "unit";
        func1 = limiters::one;
        func2 = limiters::one;
    }
    else if (contains(name, "minmod")) {
        m_name = "minmod";
        func1 = limiters::minmod;
        func2 = limiters::minmod;
    }
    else if (contains(name, "mc")) {
        m_name = "MC";
        func1 = limiters::MC;
        func2 = limiters::MC;
    }
    else if (contains(name, "leer")) {
        m_name = "van Leer";
        func1 = limiters::vanLeer;
        func2 = limiters::vanLeer;
    }
    else if (contains(name, "albada")) {
        if (contains(name, "2")) {
            m_name = "van Albada 2";
            func1 = limiters::vanAlbada2;
            func2 = limiters::vanAlbada2;
        }
        else {
            m_name = "van Albada 1";
            func1 = limiters::vanAlbada1;
            func2 = limiters::vanAlbada1;
        }
    }
    else if (contains(name, "charm")) {
        m_name = "charm";
        func1 = limiters::CHARM;
        func2 = limiters::CHARM;
    }
    else if (contains(name, "hcus")) {
        m_name = "hcus";
        func1 = limiters::HCUS;
        func2 = limiters::HCUS;
    }
    else if (contains(name, "hquick")) {
        m_name = "hquick";
        func1 = limiters::HQUICK;
        func2 = limiters::HQUICK;
    }
    else if (contains(name, "koren")) {
        m_name = "koren";
        func1 = limiters::Koren;
        func2 = limiters::Koren;
    }
    else if (contains(name, "smart")) {
        m_name = "smart";
        func1 = limiters::smart;
        func2 = limiters::smart;
    }
    else if (contains(name, "superbee")) {
        m_name = "superbee";
        func1 = limiters::superbee;
        func2 = limiters::superbee;
    }
    else {
        throw std::runtime_error("Unknown limiter '" + orig_name + "'");
    }
}

bool Limiter::is_on() const {
    return m_active;
}

const std::string& Limiter::name() const {
    return m_name;
}

double Limiter::operator()(const double& val) const {
    return func1(val);
}

double Limiter::operator()(const double& val1, const double& val2) const {
    return func2(val1, val2);
}

void Limiter::check() {
    std::set<std::string> limiters = {
            "minmod", "MC", "van Leer",
            "van Albada 1", "van Albada 2",
            "CHARM", "HCUS", "HQUICK",
            "Koren", "smart", "superbee",
    };

    double inf = std::numeric_limits<double>::infinity();
    double eps = std::numeric_limits<double>::epsilon();
    double min = std::numeric_limits<double>::min();

    std::vector<double> values = {-inf, -5.0, -1.0, -eps, -min, 0.0, min, eps, 0.5, 1.0, 2.0, 5.0, inf };

    int max_name_len = 0;
    for (const auto& name: limiters) {
        max_name_len = std::max(max_name_len, int(name.size()));
    }

    std::cout << std::setprecision(2) << std::fixed;
    std::cout << std::setw(max_name_len + 2) << " ";
    for (auto v: values) {
        std::cout << std::setw(7) << v;
    }
    std::cout << "\n";

    for (const auto& name: limiters) {
        Limiter limiter(name);

        std::cout << std::setw(max_name_len) << name << ": ";
        for (auto v: values) {
            std::cout << std::setw(7) << limiter(v);
        }
        std::cout << "\n";
    }
}

} // namespace math
} // namespace zephyr
