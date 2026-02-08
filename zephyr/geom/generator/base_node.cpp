#include <algorithm>
#include <list>
#include <utility>
#include <bits/fs_fwd.h>

#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/block.h>


namespace zephyr::geom::generator {

BaseNode::BaseNode(const Vector3d &v, bool fixed)
    : m_pos(v), m_fixed(fixed) {
}

BaseNode::Ptr BaseNode::create(const Vector3d &v, bool fixed) {
    return std::make_shared<BaseNode>(v, fixed);
}

BaseNode::Ptr BaseNode::create(double x, double y, bool fixed) {
    return std::make_shared<BaseNode>(Vector3d{x, y, 0.0}, fixed);
}

void BaseNode::set_fixed(bool fixed) {
    if (m_editable) { m_fixed = fixed; }
    std::cerr << "BaseNode::set_fixed(): BaseNode is not editable;\n";
}

void BaseNode::set_pos(const Vector3d& v) {
    if (m_editable) { m_pos = v; }
    std::cerr << "BaseNode::set_pos(): BaseNode is not editable;\n";

}

void BaseNode::set_pos(double x, double y) {
    set_pos(Vector3d{x, y, 0.0});
}

void BaseNode::reset() {
    m_editable = true;
    m_adjacent_nodes.clear();
    m_adjacent_blocks.clear();
}

// Угол блока.
// Точки упорядочены против часовой, углы тоже (+ по возрастанию)
class Angle {
public:
    Angle(Block::Ptr block_in, const BaseNode* v) : block(std::move(block_in)) {
        std::tie(v1, v2) = block->adjacent_nodes(v);
        phi1 = std::atan2(v1->y() - v->y(), v1->x() - v->x());
        phi2 = std::atan2(v2->y() - v->y(), v2->x() - v->x());
        // Углы по возрастанию
        if (phi1 > phi2) {
            phi2 += 2 * M_PI;
        }
        if (phi1 < -1.0e-3) {
            phi1 += 2 * M_PI;
            phi2 += 2 * M_PI;
        }
    }

    bool operator<(const Angle& other) const {
        return phi1 < other.phi1;
    }

    Block::Ptr block;
    BaseNode::Ptr v1, v2;
    double phi1, phi2;
};

void BaseNode::finalize(const std::set<Block::Ptr>& blocks) {
    if (blocks.empty()) {
        throw std::runtime_error("BaseNode::finalize(): No blocks were given");
    }
    if (blocks.size() == 1) {
        Block::Ptr block = *blocks.begin();
        Angle a(block, this);

        m_adjacent_nodes = {a.v1, a.v2};
        m_adjacent_blocks = {block};
        m_boundary = block->is_boundary(this);
        m_editable = false;
        return;
    }

    std::vector<Angle> angles;
    angles.reserve(blocks.size());
    for (const auto& block: blocks) {
        angles.emplace_back(block, this);
    }

    std::sort(angles.begin(), angles.end());

    std::vector<int> stops;
    for (int i = 0; i < angles.size(); ++i) {
        int j = (i - 1 + angles.size()) % angles.size();
        if (angles[j].v2 != angles[i].v1) {
            stops.emplace_back(i);
        }
    }

    // Внутренняя вершина
    if (stops.empty()) {
        for (const auto& a: angles) {
            m_adjacent_nodes.emplace_back(a.v1);
            m_adjacent_blocks.emplace_back(a.block);
            if (a.block->is_boundary(this)) {
                throw std::runtime_error("Strange Blocks: boundary condition for inner nodes");
            }
        }
        m_boundary = false;
        m_editable = false;
        return;
    }

    if (stops.size() == 1) {
        // Ровно один stop - внешняя граница
        m_adjacent_nodes.emplace_back(angles[stops.back()].v1);
        for (int k = stops.back(); k < stops.back() + angles.size(); ++k) {
            int i = k % angles.size();
            m_adjacent_nodes.emplace_back(angles[i].v2);
            m_adjacent_blocks.emplace_back(angles[i].block);
        }
        m_boundary = true;
        m_editable = false;
        return;
    }

    // stops.size() > 1. Такой случай возможен, но мне лень рассматривать
    throw std::runtime_error("BaseNode::finalize(): Strange Blocks #2");
}

bool BaseNode::regular() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::regular() called on not finalized node");
    }

    auto n_nodes = m_adjacent_nodes.size();
    auto n_blocks = m_adjacent_blocks.size();

    if (m_boundary) {
        return n_nodes == 3 && n_blocks == 2;
    }
    return n_nodes == 4 && n_blocks == 4;
}

int BaseNode::degree() const {
    return n_adjacent_nodes();
}

bool BaseNode::is_boundary() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::is_boundary() called on not finalized node");
    }
    return m_boundary;
}

const std::vector<BaseNode::WPtr>& BaseNode::adjacent_nodes() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::adjacent_nodes() called on not finalized node");
    }
    return m_adjacent_nodes;
}

int BaseNode::n_adjacent_nodes() const {
    return static_cast<int>(adjacent_nodes().size());
}

const std::vector<Block::WPtr> &BaseNode::adjacent_blocks() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::adjacent_blocks() called on not finalized node");
    }
    return m_adjacent_blocks;
}

int BaseNode::n_adjacent_blocks() const {
    return static_cast<int>(adjacent_blocks().size());
}

} // namespace zephyr::geom::generator
