#pragma once
#include <vector>
#include <limits>
#include <stdexcept>

#include <zephyr/utils/json.h>

namespace zephyr::io {

class WriteFreq {
public:

    explicit WriteFreq(const utils::Json& config, double init_time = 0.0);

    bool time_to_write(double timestep);

protected:
    std::vector<double> m_timecodes;
    std::vector<double> m_frequencies;
    double m_next_write;
    double m_time_delta;
};

} // namespace zephyr::io