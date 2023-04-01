#include <iomanip>

#include <zephyr/utils/stopwatch.h>

namespace zephyr { namespace utils {

Stopwatch::Stopwatch(bool run)
        : m_up(false), m_start(), m_elapsed(duration::zero()) {
    if (run) {
        start();
    }
}

Stopwatch::duration Stopwatch::elapsed() const {
    if (!m_up) {
        return m_elapsed;
    }
    else {
        return m_elapsed + (clock::now() - m_start);
    }
}

void Stopwatch::start() {
    m_up = true;
    m_elapsed = duration::zero();
    m_start = clock::now();
}

void Stopwatch::stop() {
    if (m_up) {
        m_elapsed += (clock::now() - m_start);
        m_up = false;
    }
}

void Stopwatch::resume() {
    if (!m_up) {
        m_start = clock::now();
        m_up = true;
    }
}

bool Stopwatch::is_up() const {
    return m_up;
}

long Stopwatch::milliseconds() const {
    return int(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed()).count());
}

long Stopwatch::seconds() const {
    return int(std::chrono::duration_cast<std::chrono::seconds>(elapsed()).count());
}

long Stopwatch::minutes() const {
    return int(std::chrono::duration_cast<std::chrono::minutes>(elapsed()).count());
}

long Stopwatch::hours() const {
    return std::chrono::duration_cast<std::chrono::hours>(elapsed()).count();
}

long Stopwatch::days() const {
    return hours() / 24;
}

std::string Stopwatch::extended_time() const {
    long full = seconds();

    long minutes = full / 60 % 60;
    long hours = full / 3600 % 24;
    long days = full / 86400;
    long seconds = full % 60;

    std::stringstream ss;
    if (days > 0) {
        ss << std::setw(2) << days << " d ";
    } else {
        ss << "     ";
    }
    if (hours > 0 || days > 0) {
        ss << std::setw(2) << hours << " h ";
    } else {
        ss << "     ";
    }
    ss << std::setw(2) << minutes << " m ";
    ss << std::setw(2) << seconds << " s";

    return ss.str();
}

} // utils
} // zephyr