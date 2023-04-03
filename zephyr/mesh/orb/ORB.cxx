#include <zephyr/network/decomposition/ORB.h>

namespace zephyr { namespace network { namespace decomposition {

inline double fit(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

ORB::ORB(
    Network&  network,
    Storage&  elements,
    Domain&   domain,
    Measurer& measurer,
    if_multithreading(ThreadPool &pool,)
    bool      newton,
    double    mobility
)
    : Base(network, elements, domain, measurer if_multithreading(, pool)),
      m_newton(newton),
      m_mobility(fit(mobility, 0.0, 0.45))
{


}

#ifdef ZEPHYR_ENABLE_YAML

ORB::ORB(
    Network&  network,
    Storage&  elements,
    Domain&   domain,
    Measurer& measurer,
    const YAML::Node &config
    if_multithreading(, ThreadPool &pool)
)
    : Base(network, elements, domain, measurer if_multithreading(, pool)),
      m_newton(defaults::newton),
      m_mobility(defaults::mobility)
{

    if (config["newton"]) {
        use_newton(config["newton"].as<bool>());
    }

    if (config["mobility"]) {
        set_mobility(config["mobility"].as<double>());
    }

    if (!config["blocks"]) {
        std::string message = "Декомпозиция ORB. "
                              "Необходимо указать ключ 'blocks' внутри 'decomposition' "
                              "с порядком декомпозиции по координатам. Возможные значения: "
                              "'X', 'Y', 'Z', 'XY', 'YX', 'R', 'P', 'RP' и так далее";
        if (network.is_master()) {
            std::cerr << message << "\n";
        }
        throw std::runtime_error(message);
    }

    std::string type = config["blocks"].as<std::string>();

    int n_proc = network.size();
    int n_proc_1 = config["n_proc_1"] ? config["n_proc_1"].as<int>() : -1;
    int n_proc_2 = config["n_proc_2"] ? config["n_proc_2"].as<int>() : -1;
    int n_proc_3 = config["n_proc_3"] ? config["n_proc_3"].as<int>() : -1;

    if (n_proc_1 < 0 && n_proc_2 < 0 && n_proc_3 < 0) {
        // Вариант 1. Автоматический конструктор по геометрическим размерам.
        // Внимание: не работает с полярным разбиением.

        Vector3d sizes = domain.get_bounding_box_sizes();
        blocks = std::unique_ptr<Blocks>(new Blocks(type, n_proc, sizes));
    } else {
        // Вариант 2. С полным или частичным заданием размеров.
        // Для двумерного разбиения требуется указать размеров хотя бы
        // по одной из осей
        // Для трехмерного разбиения требуется указать размеры хотя бы
        // по двум осям
        // Для одномерного разбиения указание размеров не требуется.

        blocks = std::unique_ptr<Blocks>(new Blocks(type, n_proc, n_proc_1, n_proc_2, n_proc_3));
    }

    Vector3d dc = domain.get_bounding_box_center();
    Vector3d ds = domain.get_bounding_box_sizes();

    blocks->setup_bounds(dc - 0.5 * ds, dc + 0.5 * ds);

    //if (network.is_master()) {
    //    blocks->info();
    //}

    r_search.resize(blocks->size(), 0.0);

    redistribute();
}

#endif

void ORB::use_newton(bool val) {
    m_newton = val;
}

void ORB::set_mobility(double val) {
    m_mobility = fit(val, 0.0, 0.45);
}

double ORB::mobility() const {
    return m_mobility;
}

bool ORB::newton() const {
    return m_newton;
}

Vector3d ORB::center() {
    return blocks->center(m_rank);
}

void ORB::balancing(const std::vector<double> &w) {
    if (m_newton) {
        blocks->balancing_newton(w, m_mobility);
    }
    else {
        blocks->balancing_simple(w, m_mobility);
    }
}

int ORB::rank(const Vector3d& v) const {
    return blocks->rank(v);
}

bool ORB::is_near(const Vector3d& v, int neib_rank) const {
    double search_radius = std::max(
            r_search[m_rank],
            r_search[neib_rank]
    );

    return blocks->is_near(m_rank, neib_rank, v, search_radius);
}

void ORB::collect_local_info() {
    using data::neibsSearchRadius;

    auto func = [](Storage::Item elem) -> double {
        return elem[neibsSearchRadius].value;
    };

    double nR = 0.0;
#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
       nR = m_pool.max(m_locals.begin(), m_locals.end(), 0.0, func);
    }
    else
#endif
    {
        for (auto elem: m_locals) {
            nR = std::max(nR, func(elem));
        }
    }

    net.all_gather(nR, r_search);
}

void ORB::collect_aliens_info() {
    // Класс не содержит никаких данных, которые строятся
    // по спискам aliens.
}

} // decomposition
} // network
} // zephyr