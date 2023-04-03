#include <fstream>
#include <numeric>
#include <algorithm>

#include <zephyr/network/decomposition/base.h>

#include <zephyr/network/decomposition/VD3.h>
#include <zephyr/network/decomposition/ORB.h>

namespace zephyr { namespace network { namespace decomposition {

Base::Base(
	Network&  net,
	Storage&  elements,
	Domain&   domain,
	Measurer& measurer
	if_multithreading(, ThreadPool& pool)
)
    : net(net),
      if_multithreading(m_pool(pool),)
      m_locals(elements),
      m_domain(domain),
      m_measurer(&measurer),
      m_rank(net.rank())
{
    init_neibs();
    init_tourists();
    init_migrants();

    has_exchange_layer = m_locals.has_type(neibsSearchRadius);
}

#ifdef ZEPHYR_ENABLE_YAML
std::unique_ptr<Base> Base::create(
    Network&  net,
    Storage&  elements,
    Domain&   domain,
    Measurer& measurer,
    if_multithreading(ThreadPool& pool,)
    const YAML::Node& config
) {
    if (!config["type"]) {
        throw std::runtime_error("Please setup decomposition type");
    }
    std::string type = config["type"].as<std::string>();

    if (type == "VD3") {
        return std::unique_ptr<decomposition::VD3>(
                new decomposition::VD3(net, elements, domain, measurer, config if_multithreading(, pool)));
    }
    else if (type == "ORB") {
        return std::unique_ptr<decomposition::ORB>(
                new decomposition::ORB(net, elements, domain, measurer, config if_multithreading(, pool)));
    } else {
        throw std::runtime_error("Unknown decomposition type '" + type + "'");
    }
}
#endif


void Base::init_neibs() {
    m_neibs.resize(net.size());
    for (int i = 0; i < int(net.size()); ++i) {
        m_neibs[i] = i;
    }
}

void Base::init_migrants() {
	migrants_out.clear(); migrants_out.resize(net.size());
	migrants_in.clear();  migrants_in.resize(net.size());

	for (auto& i: migrants_out) i.init(m_locals.type_names(), {0});
	for (auto& i: migrants_in)  i.init(m_locals.type_names(), {0});
}

void Base::init_tourists() {
	m_tourists_labels.clear();

	auto size = net.size();
	m_tourists_labels.resize(size);

	for (auto i = 0; i < size; i++) {
		m_tourists_labels[i].clear();
	}

	m_tourists.init(m_locals.type_names(), {0});
	m_aliens.init(m_locals.type_names(), {0});
}

enum Tags {
    DATA_TAG
};

void Base::exchange_start() {
    if (!has_exchange_layer) {
        return;
    }

    m_send_timer.resume();

    m_send1_timer.resume();

    size_t counter = 0;
    {
        using namespace zephyr::data;
        for (std::size_t l = 0; l < m_tourists_labels.size(); ++l) {
            for (auto label: m_tourists_labels[l]) {
                m_tourists[counter][item] = m_locals[label][item];
                ++counter;
            }
        }
    }
    m_send1_timer.stop();

    m_send2_timer.resume();

    m_send_request.resize(net.size());
    m_recv_request.resize(net.size());

    size_t offset1 = 0;
    for (int rank = 0; rank < int(net.size()); ++rank) {
        if (m_send_count[rank] > 0) {
            char* data = m_tourists.data() + offset1 * m_tourists.itemsize();
            MPI_Isend(data, m_send_count[rank] * m_tourists.itemsize(), MPI_CHAR, rank,
                      DATA_TAG, net.comm(), &m_send_request[rank]);
            offset1 += m_send_count[rank];
        }
    }

    size_t offset2 = 0;
    for (int rank = 0; rank < int(net.size()); ++rank) {
        if (m_recv_count[rank] > 0) {
            char* data = m_aliens.data() + offset2 * m_aliens.itemsize();
            MPI_Irecv(data, m_recv_count[rank] * m_aliens.itemsize(),
                      MPI_CHAR, rank, DATA_TAG, net.comm(), &m_recv_request[rank]);

            offset2 += m_recv_count[rank];
        }
    }

    m_send2_timer.stop();

    m_send_timer.stop();
}

void Base::exchange_end() {
    if (!has_exchange_layer) {
        return;
    }

    m_wait_timer.resume();

    m_measurer->stop();

    m_send_status.resize(net.size());
    for (int rank = 0; rank < int(net.size()); ++rank) {
        if (m_send_count[rank] > 0) {
            MPI_Wait(&m_send_request[rank], &m_send_status[rank]);
        }
    }

    m_recv_status.resize(net.size());
    for (int rank = 0; rank < int(net.size()); ++rank) {
        if (m_recv_count[rank] > 0) {
            MPI_Wait(&m_recv_request[rank], &m_recv_status[rank]);
        }
    }

    m_wait_timer.stop();
    m_recv_timer.resume();

    m_measurer->resume();

    remove_period_for(m_aliens);

    m_recv_timer.stop();
}

void Base::exchange() {
    exchange_start();
    exchange_end();
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
static std::mutex mutex1;
#endif

void Base::mark_tourists_out(Storage::Item it, Storage::Item begin) {
    for (int neib_rank = 0; neib_rank < int(net.size()); ++neib_rank) {
        if (neib_rank == m_rank) {
            continue;
        }

        if (is_near((Vector3d &) it[coords], neib_rank)) {
#ifndef ZEPHYR_ENABLE_MULTITHREADING
            m_tourists_labels[neib_rank].push_back(it - begin);
#else
            mutex1.lock();
            m_tourists_labels[neib_rank].push_back(it - begin);
            mutex1.unlock();
#endif
        }
    }
}

void Base::prepare_aliens() {
    if (!has_exchange_layer) {
        return;
    }

    m_prep_timer.resume();

    for (auto& i: m_tourists_labels) i.clear();

    auto b = m_locals.begin();

    auto func = [&](Storage::Item it) {
        mark_tourists_out(it, b);
    };

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
        m_pool.for_each(m_locals.begin(), m_locals.end(), func);
    }
    else
#endif
    {
        for (auto& it : m_locals) {
            func(it);
        }
    }

    // Обновить список соседних процессов
    m_neibs.clear();
    for (int rank = 0; rank < net.size(); rank++) {
        if (!m_tourists_labels[rank].empty()) {
            m_neibs.push_back(rank);
        }
    }

    m_send_count.resize(net.size());
    for (int r = 0; r < int(net.size()); ++r) {
        m_send_count[r] = m_tourists_labels[r].size();
    }

    net.all_to_all(m_send_count, m_recv_count);

    size_t n_aliens = std::accumulate(m_recv_count.begin(), m_recv_count.end(), 0ul);

    m_aliens.resize(n_aliens);

    size_t n_tourists = std::accumulate(m_send_count.begin(), m_send_count.end(), 0ul);

    m_tourists.resize(n_tourists);

    m_prep_timer.stop();
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
static std::mutex mutexb1;
static std::mutex mutexb2;
#endif
void Base::prepare_migrants() {
    auto b = m_locals.begin();
#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
        auto func = [&](Storage::Item it) {
            auto s = rank((Vector3d&) it[coords]);
            if (s != m_rank) {
                mutexb1.lock(); migrants_out[s].append(it[item]); mutexb1.unlock();
                mutexb2.lock(); to_remove.emplace_back(it - b); mutexb2.unlock();
            }
        };
        m_pool.for_each(m_locals.begin(), m_locals.end(), func);
    }
    else
#endif
    {
        for (auto it : m_locals) {
            auto s = rank((Vector3d &) it[coords]);
            if (s != m_rank) {
                migrants_out[s].append(it[item]);
                to_remove.emplace_back(it - b);
            }
        }
    }
}

void Base::exchange_migrants() {
    for (auto& i: migrants_out) i.resize(0);
    for (auto& i: migrants_in)  i.resize(0);

    to_remove.resize(0);

    fit_in_period_for(m_locals);
    prepare_migrants();

    for (std::size_t i = 0; i < net.node.size(); i++) {
        for (std::size_t j = 0; j < net.node.size(); j++) {
            net.node[i].async(migrants_in[j]) =
                    net.node[j].async(migrants_out[i]);
        }
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
        std::sort(to_remove.begin(), to_remove.end());
    }
#endif

    m_locals.remove(to_remove);

    for (std::size_t i = 0; i < net.node.size(); i++)
        net.node[i].await();

    for (auto& migrants : migrants_in)
        m_locals.join(migrants);

    remove_period_for(m_locals);
}

void Base::redistribute() {
    exchange_migrants();
    collect_local_info();
    prepare_aliens();
    exchange();
    collect_aliens_info();
}

inline load_t imbalance(const std::vector<load_t>& loads) {
    load_t sum = 0.0;
    load_t max = -1.0;
    for (auto w: loads) {
        sum += w;
        max = std::max(w, max);
    }
    load_t avg = sum / loads.size();

    return max / avg - 1.0;
}

void Base::update() {
    auto w = m_measurer->get_all(net);
    if (w.size() != net.size()) {
        std::string message = "void Base::update() error: sizes of storage and w are not equal";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    if (net.is_master()) {
        static size_t counter = 0;
        static std::ofstream file;
        static auto start = std::chrono::steady_clock::now();

        file.open("imbalance.txt", counter < 1 ? std::ios_base::out : std::ios_base::app);

        auto finish = std::chrono::steady_clock::now();
        auto seconds = long(std::chrono::duration<double>(finish - start).count());

        file << counter << " " << seconds << " " << imbalance(w) << "\n";

        file.close();

        ++counter;
    }

    balancing(w);
    redistribute();
    m_measurer->start();
}

void Base::permute_tourists_labels_out(std::vector<std::size_t>& perm) {
	for (auto& list : m_tourists_labels) {
        for (auto &i : list) {
            i = perm[i];
        }
    }
}

void Base::collect_at(std::size_t rank) {
	for (auto& i: migrants_out) i.resize(0);
	for (auto& i: migrants_in)  i.resize(0);
	
	if (m_rank != rank) {
		migrants_out[rank] = m_locals;
		m_locals.resize(0);
	}

	for (std::size_t i = 0; i < net.node.size(); ++i) {
		net.node[rank](migrants_in[i]) = 
			net.node[i](migrants_out[rank]);
	}

	for (auto& migrants : migrants_in) {
		m_locals.join(migrants);
	}
}

Storage& Base::inner_elements() { return m_locals; }
Storage& Base::outer_elements() { return m_aliens; }

void Base::signal_to_permute(std::vector<std::size_t>& perm) {
    permute_tourists_labels_out(perm);
};

void Base::set_load_measurer(Measurer& lm) {
    m_measurer = &lm;
    m_measurer->start();
}

Measurer& Base::load_measurer() const {
    return *m_measurer;
}

Network& Base::network() const { return net; }

Storage Base::average(
	std::vector<std::string> fields,
	std::pair<long long int*, std::size_t> indices
) {

	using ::zephyr::data::density;
	using ::zephyr::data::pressure;
	using ::zephyr::data::uid;

	auto selected = m_locals[fields];
	if (indices.first != nullptr) {
		selected = selected[indices];
	}
	auto size = selected.size();

	std::vector<Storage> averages(net.size());
	fields.push_back("uid");
	averages[m_rank] = Storage(fields, 1);
	auto& avg = averages[m_rank];
	avg[0][uid].value = selected.size();
	auto itemSize = selected.itemsize();
	auto scalarCount = itemSize/sizeof(double);

	auto aData = avg.data();
	auto eData = selected.data();

	for (std::size_t i = 0; i < size; ++i) {
		for (std::size_t scalarIndex = 0; scalarIndex < scalarCount; ++scalarIndex) {
			double& aScalar = *(((double*) aData) + scalarIndex);
			double& eScalar = *(((double*) (eData + i*itemSize)) + scalarIndex);
			aScalar += eScalar;
		}
	}

	for (std::size_t i = 0; i < net.node.size(); i++) {
		for (std::size_t j = 0; j < net.node.size(); j++) {
			if (i != j)
				net.node[i](averages[j]) = 
					net.node[j](averages[j]);
		}
	}

	for (std::size_t scalarIndex = 0; scalarIndex < scalarCount; ++scalarIndex) {
		double total = 0.0;
		std::size_t count = 0;
		for (std::size_t i = 0; i < net.node.size(); i++) {
			auto aData = averages[i].data();
			double& aScalar = *(((double*) aData) + scalarIndex);
			total += aScalar;
			count += averages[i][0][uid];
		}

		auto aData = averages[m_rank].data();
		double& aScalar = *(((double*) aData) + scalarIndex);
		aScalar = total/count;
	}

	fields.pop_back();
	return averages[m_rank][fields];
}

Storage Base::average(
	std::vector<std::string> fields,
	std::vector<long long int>& index
) {
	return average(fields, {index.data(), index.size()});
}

void Base::remove_period(Storage::Item element) {
    auto domain_center = center();
    auto& element_coords = (Vector3d&) element[coords];
    auto offset = m_domain.shortest(element_coords - domain_center);
    element_coords = domain_center + offset;
}

void Base::fit_in_period(Storage::Item element) {
    auto& element_coords = (Vector3d&) element[coords];
    m_domain.fit_in_period(element_coords);
}

void Base::fit_in_period_for(Storage& s) {
    auto func = [&](Storage::Item it) {
        fit_in_period(it);
    };
#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
        m_pool.for_each(s.begin(), s.end(), func);
    }
    else
#endif
    {
        for (auto &it : s) {
            func(it);
        }
    }
}

void Base::remove_period_for(Storage& s) {
    auto func = [&](Storage::Item it) {
        remove_period(it);
    };
#ifdef ZEPHYR_ENABLE_MULTITHREADING
    if (m_pool.is_active()) {
        m_pool.for_each(s.begin(), s.end(), func);
    }
    else
#endif
    {
        for (auto &it : s) {
            func(it);
        }
    }
}

void Base::to_global_CS() {
    fit_in_period_for(inner_elements());
    fit_in_period_for(outer_elements());
}

void Base::to_local_CS() {
    remove_period_for(inner_elements());
    remove_period_for(outer_elements());
}


} /*decomposition*/ } /*network*/ } /*zephyr*/