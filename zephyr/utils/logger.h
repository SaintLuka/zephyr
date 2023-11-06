#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <fstream>
#include <boost/variant.hpp>
#include <zephyr/math/cfd/models.h>

namespace zephyr::utils {

using namespace std;

struct Output : public boost::static_visitor<> {
    std::ostream &os;

    Output(std::ostream &os) : os(os) {}

    template<typename T>
    void operator()(const T &v) const {
        os << v;
    }

    ~Output() = default;
};

class Logger {
public:
    using types = boost::variant<int, double, size_t, long long,
            math::mmf::PState, math::mmf::QState, math::mmf::Flux,
            math::smf::PState, math::smf::QState, math::smf::Flux,
            phys::Fractions, phys::FractionsFlux>;

    Logger();

    void log(const string &name, const types &value);

    void save_to_file(const string &filename = "../../problems/logs.csv");

    vector<pair<string, types>> logs;
};

}
