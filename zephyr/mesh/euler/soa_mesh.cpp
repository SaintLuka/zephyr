#include <fstream>

#include <zephyr/mesh/euler/soa_mesh.h>

#include <zephyr/mesh/amr2/apply.h>
#include <zephyr/mesh/amr2/balancing.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/mesh/euler/soa_prim.h>

using namespace zephyr::geom;
using generator::Rectangle;
using generator::Cuboid;

namespace zephyr::mesh {




SoaMesh::SoaMesh(Generator& gen) {
    gen.initialize(m_locals);

    structured = false;
    m_nx = 1;
    m_ny = 1;
    m_nz = 1;

    if (dynamic_cast<Rectangle*>(&gen) != nullptr) {
        Rectangle* rect = dynamic_cast<Rectangle*>(&gen);
        if (rect->structured()) {
            structured = true;
            m_nx = rect->nx();
            m_ny = rect->ny();
            m_nz = 1;
        }
    }
    else if (dynamic_cast<Cuboid*>(&gen) != nullptr) {
        Cuboid* cuboid = dynamic_cast<Cuboid*>(&gen);
        structured = true;
        m_nx = cuboid->nx();
        m_ny = cuboid->ny();
        m_nz = cuboid->nz();
    }
}

void SoaMesh::init_amr() {
    // Эта функция нифига не используется же?
    m_max_level = 0;

    if (m_locals.empty()) {
        return;
    }

    int rank = mpi::rank();
    int size = mpi::size();

    auto cells_nums = mpi::all_gather(m_locals.size());
    std::vector<int> offset(size, 0);
    for (int r = 0; r < size - 1; ++r) {
        offset[r + 1] = offset[r] + cells_nums[r];
    }

    for (int ic = 0; ic < m_locals.size(); ++ic) {
        m_locals.b_idx[ic] = offset[rank] + ic;
        m_locals.z_idx[ic] = 0;
        m_locals.level[ic] = 0;
        m_locals.flag[ic] = 0;
    }
}

bool SoaMesh::is_adaptive() const {
    return m_max_level > 0;
}

int SoaMesh::max_level() const {
    return m_max_level;
}

void SoaMesh::set_max_level(int max_level) {
    m_max_level = std::max(0, std::min(max_level, 15));
}

void SoaMesh::set_distributor(const std::string& name) {
    if (name == "empty") {
        distributor = Distributor::empty();
    } else {
        distributor = Distributor::simple();
    }
}

void SoaMesh::set_distributor(Distributor distr) {
    distributor = std::move(distr);
}

void SoaMesh::balance_flags() {
#if SCRUTINY
    static bool first_time = true;
    if (first_time) {
        if (check_base() < 0) {
            throw std::runtime_error("Check base failed");
        }
        first_time = false;
    }
#endif
    if constexpr (mpi::single()) {
        amr2::balance_flags(m_locals, m_max_level);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::balance_flags(m_locals, m_aliens, m_max_level, *this);
    }
#endif
}

void SoaMesh::apply_flags() {
    if constexpr (mpi::single()) {
        amr2::apply(m_locals, distributor);
    }
#ifdef ZEPHYR_MPI
    else {
        amr2::apply(m_locals, m_aliens, distributor, *this);
    }
#endif

#if SCRUTINY
    if (check_refined() < 0) {
        //throw std::runtime_error("Check refined failed");
    }
#endif
}

void SoaMesh::refine() {
    if (!is_adaptive()) { return; }

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::single() && m_locals.empty()) {
        throw std::runtime_error("EuMesh::refine() error: Empty mesh");
    }

    full.resume();

    balance.resume();
    balance_flags();
    balance.stop();

    apply.resume();
    apply_flags();
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr2::check_frequency == 0) {
        std::cout << "  Balance steps elapsed: " << std::setw(10) << balance.milliseconds() << " ms\n";
        std::cout << "  Apply steps elapsed:   " << std::setw(10) << apply.milliseconds() << " ms\n";
        std::cout << "Refine steps elapsed:    " << std::setw(10) << full.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}


int SoaMesh::check_base() const {
    if (m_locals.empty()) {
        if constexpr (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals.dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (index_t ic = 0; ic < m_locals.size(); ++ic) {
        if (m_locals.index[ic] < 0 || m_locals.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (m_locals.rank[ic] < 0 || m_locals.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (m_locals.faces.is_undefined(m_locals.face_begin[ic] + i)) {
                std::cout << "\tCell has no one of main faces\n";
                m_locals.print_info(ic);
                return -1;
            }
        }
        if (m_locals.face_count(ic) > FpC(dim)) {
            std::cout << "\tCell has too much faces (" << m_locals.face_count(ic) << ")\n";
            m_locals.print_info(ic);
            return -1;
        }

        // Число вершин ???

        // Правильное задание геометрии
        res = m_locals.check_geometry(ic);
        if (res < 0) return res;

        // Грани правильно ориентированы
        res = m_locals.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = m_locals.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка смежности
        res = m_locals.check_connectivity(ic);
        if (res < 0) return res;
    }

    return 0;
}

int SoaMesh::check_refined() const {
    if (m_locals.empty()) {
        if constexpr (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals.dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (index_t ic = 0; ic < m_locals.size(); ++ic) {
        if (m_locals.is_undefined(ic)) {
            std::cout << "\tUndefined cell\n";
            return -1;
        }

        if (m_locals.index[ic] < 0 || m_locals.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (m_locals.rank[ic] < 0 || m_locals.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (m_locals.faces.is_undefined(m_locals.face_begin[ic] + i)) {
                std::cout << "\tCell has no one of main faces\n";
                m_locals.print_info(ic);
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = m_locals.node_begin[ic]; i < m_locals.node_begin[ic + 1]; ++i) {
            for (int j = i + 1; j < m_locals.node_begin[ic + 1]; ++j) {
                double dist = (m_locals.verts[i] - m_locals.verts[j]).norm();
                if (dist < 1.0e-5 * m_locals.linear_size(ic)) {
                    std::cout << "\tIdentical vertices\n";
                    m_locals.print_info(ic);
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = m_locals.check_geometry(ic);
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = m_locals.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = m_locals.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка сложных граней
        res = m_locals.check_complex_faces(ic);
        if (res < 0) return res;

        // Проверка смежности
        res = m_locals.check_connectivity(ic);
        if (res < 0) return res;
    }

    return 0;
}

Box SoaMesh::bbox() const {
    Box box1 = Box::Empty(3);
    for (auto& v: m_locals.verts) {
        box1.capture(v);
    }

    Box box2(box1);

#ifdef ZEPHYR_MPI
    if (!mpi::single()) {
        // Покомпонентный минимум/максимум
        MPI_Allreduce(box1.vmin.data(), box2.vmin.data(), 3, MPI_DOUBLE, MPI_MIN, mpi::comm());
        MPI_Allreduce(box1.vmax.data(), box2.vmax.data(), 3, MPI_DOUBLE, MPI_MAX, mpi::comm());
    }
#endif

    return box2;
}

void SoaMesh::push_back(const geom::Line& line) {
    geom::Polygon poly = {line[0], line[1], line[1], line[0]};
    m_locals.push_back(poly);
}

void SoaMesh::push_back(const geom::Polygon& poly) {
    m_locals.push_back(poly);
}

void SoaMesh::push_back(const geom::Polyhedron& poly) {
    throw std::runtime_error("Polyhedrpushback");
}

CellIt SoaMesh::begin() { return {&m_locals, 0, &m_aliens}; }

CellIt SoaMesh::end() { return {&m_locals, m_locals.size(), &m_aliens}; }


QCell SoaMesh::operator()(int i, int j) {
    i = (i + m_nx) % m_nx;
    j = (j + m_ny) % m_ny;
    return operator[](m_ny * i + j);
}

QCell SoaMesh::operator()(int i, int j, int k) {
    i = (i + m_nx) % m_nx;
    j = (j + m_ny) % m_ny;
    k = (k + m_nz) % m_nz;
    return operator[](m_nz * (m_ny * i + j) + k);
}

} // namespace zephyr::mesh