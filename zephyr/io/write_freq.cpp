#include <algorithm>
#include <zephyr/io/write_freq.h>

namespace zephyr::io {

inline constexpr double inf() {
    return std::numeric_limits<double>::infinity();
}

WriteFreq::WriteFreq(const utils::Json &config, double init_time) {
    auto add_one = [this](const utils::Json& config) {
        if (config["time"]) {
            m_timecodes.push_back(config["time"].as<double>());
        } else {
            throw std::runtime_error("WriteFreq::WriteFreq: not found 'time'");
        }

        if (config["freq"]) {
            m_frequencies.push_back(config["freq"].as<double>());
        } else if (config["frequency"]) {
            m_frequencies.push_back(config["frequency"].as<double>());
        } else {
            throw std::runtime_error("WriteFreq::WriteFreq: not found 'freq' or 'frequency'");
        }
    };

    if (config.is_object()) {
        add_one(config);
    }
    else if (config.is_array()) {
        for (auto& elem: config.array_items()) {
            add_one(elem);
        }
    }
    else {
        throw std::runtime_error("Can't initialize write_dreq");
    }

    double min_freq = *std::min_element(m_frequencies.begin(), m_frequencies.end());

    m_next_write = min_freq * static_cast<unsigned int>(init_time / min_freq);
    m_time_delta = std::numeric_limits<double>::max();
}

bool WriteFreq::time_to_write(double timestep) {
    bool res = timestep >= m_next_write;

    if (res) {
        for (size_t i = 0; i < m_timecodes.size(); ++i) {
            if (timestep >= m_timecodes[i]) {
                double next_timecode = i < m_timecodes.size() - 1 ? m_timecodes[i + 1] : inf();
                m_time_delta = std::min(m_frequencies[i], next_timecode - m_next_write);
            }
        }
        m_next_write += m_time_delta;
    }
    return res;
}

} // namespace zephyr::io