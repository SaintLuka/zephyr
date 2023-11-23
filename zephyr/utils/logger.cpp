#include <zephyr/utils/logger.h>

namespace zephyr::utils {

void Logger::log(const string &name, const types &value) {
    logs.emplace_back(name, value);
}

void Logger::save_to_file(const string &filename) {
    std::ofstream out(filename, std::ios::trunc);
    if (!out.is_open())
        throw runtime_error("can't open file for logs");

    Output output(out);

    for (size_t i = 0; i < logs.size(); i++) {
        out << logs[i].first << "; ";
        boost::apply_visitor(output, logs[i].second);
        out << "\n";
    }

    out.close();
}

Logger::Logger() : logs() {}

}