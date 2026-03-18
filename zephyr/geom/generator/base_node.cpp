#include <algorithm>

#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/block.h>

namespace zephyr::geom::generator {

BaseNode::BaseNode(const Vector3d &v, bool fixed)
    : m_fixed(fixed), m_pos(v) {
}

BaseNode::Ptr BaseNode::create(const Vector3d &v, bool fixed) {
    return std::make_shared<BaseNode>(v, fixed);
}

BaseNode::Ptr BaseNode::create(double x, double y, bool fixed) {
    return std::make_shared<BaseNode>(Vector3d{x, y, 0.0}, fixed);
}

void BaseNode::set_fixed(bool fixed) {
    if (m_editable) { m_fixed = fixed; }
    else {
        std::cerr << "BaseNode::set_fixed warning: BaseNode is not editable;\n";
    }
}

void BaseNode::set_pos(double x, double y) {
    set_pos(Vector3d{x, y, 0.0});
}

void BaseNode::set_pos(const Vector3d& v) {
    if (m_editable) { m_pos = v; }
    else {
        std::cerr << "BaseNode::set_pos warning: BaseNode is not editable;\n";
    }
}

void BaseNode::clear() {
    m_editable = true;
    m_adjacent_nodes.clear();
    m_adjacent_blocks.clear();
}

// Угол блока. Точки упорядочены против часовой стрелки, углы положительные
// и упорядочены по возрастанию (угол поворота относительно основной вершины).
class Corner {
public:
    Corner(Block::Ptr block_in, const BaseNode* v)
        : block(std::move(block_in)) {
        // adjacent_nodes гарантирует порядок обхода
        std::tie(v1, v2) = block->adjacent_nodes(v);
        phi1 = std::atan2(v1->y() - v->y(), v1->x() - v->x());
        phi2 = std::atan2(v2->y() - v->y(), v2->x() - v->x());
        // Расположим углы по возрастанию
        if (phi2 < phi1) {
            phi2 += 2 * M_PI;
        }
        // Сделаем оба угла положительными
        if (phi1 < -1.0e-3) {
            phi1 += 2 * M_PI;
            phi2 += 2 * M_PI;
        }
    }

    // Позволяет отсортировать блоки против часовой стрелки
    bool operator<(const Corner& other) const {
        return phi1 < other.phi1;
    }

    Block::Ptr block;
    BaseNode::Ptr v1, v2;
    double phi1, phi2;
};

void BaseNode::finalize(const std::set<Block::Ptr>& blocks) {
    if (!m_editable) {
        throw std::runtime_error("BaseNode::finalize: BaseNode is not editable");
    }
    if (blocks.empty()) {
        throw std::runtime_error("BaseNode::finalize: no blocks were given");
    }
    if (blocks.size() == 1) {
        Block::Ptr block = *blocks.begin();
        auto[v1, v2] = block->adjacent_nodes(this);
        m_adjacent_nodes = {v1, v2};
        m_adjacent_blocks = {block};
        if (auto[side1, side2] = block->incident_sides(this);
            !block->boundary(side1) || !block->boundary(side2)) {
            throw std::runtime_error("BaseNode::finalize: invalid corner BaseNode");
        }
        m_boundary = true;
        m_editable = false;
        return;
    }

    // Собрать уголки блоков
    std::vector<Corner> corners;
    corners.reserve(blocks.size());
    for (const auto& block: blocks) {
        corners.emplace_back(block, this);
    }

    // Сортирует блоки против часовой стрелки
    std::sort(corners.begin(), corners.end());

    // "Разрывы" в порядке блоков, достигаются на границе области
    std::vector<int> stops;
    for (int i = 0; i < corners.size(); ++i) {
        int j = (i - 1 + corners.size()) % corners.size();
        if (corners[j].v2 != corners[i].v1) {
            stops.emplace_back(i);
        }
    }

    // Разрывов нет, это внутренняя вершина
    if (stops.empty()) {
        for (const auto& a: corners) {
            m_adjacent_nodes.emplace_back(a.v1);
            m_adjacent_blocks.emplace_back(a.block);

            auto [side1, side2] = a.block->incident_sides(this);
            if (a.block->boundary(side1) || a.block->boundary(side2)) {
                throw std::runtime_error("BaseNode::finalize: strange blocks, boundary condition for inner nodes");
            }
        }
        m_boundary = false;
        m_editable = false;
        return;
    }

    // Один разрыв, это нормальная граничная вершина
    if (stops.size() == 1) {
        // Ровно один stop - внешняя граница
        m_adjacent_nodes.emplace_back(corners[stops.back()].v1);
        for (int k = stops.back(); k < stops.back() + corners.size(); ++k) {
            int i = k % corners.size();
            m_adjacent_nodes.emplace_back(corners[i].v2);
            m_adjacent_blocks.emplace_back(corners[i].block);
        }
        m_boundary = true;
        m_editable = false;
        return;
    }

    // stops.size() > 1. Такой случай возможен, но мне лень рассматривать
    // Пример: у вершины два смежных блока, которые друг для друга смежными
    // не являются. Получается нетрадиционная граница.
    throw std::runtime_error("BaseNode::finalize: strange blocks, blocks is adjacent only through BaseNode");
}

bool BaseNode::boundary() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::is_boundary: called for not finalized BaseNode");
    }
    return m_boundary;
}

bool BaseNode::regular() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::regular: called for not finalized BaseNode");
    }
    auto n_nodes  = n_adjacent_nodes();
    auto n_blocks = n_adjacent_blocks();
    if (m_boundary) {
        return n_nodes == 3 && n_blocks == 2;
    }
    return n_nodes == 4 && n_blocks == 4;
}

int BaseNode::degree() const {
    return n_adjacent_nodes();
}

int BaseNode::n_adjacent_nodes() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::n_adjacent_nodes: called for not finalized BaseNode");
    }
    return static_cast<int>(adjacent_nodes().size());
}

int BaseNode::n_adjacent_blocks() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::n_adjacent_blocks: called for not finalized BaseNode");
    }
    return static_cast<int>(adjacent_blocks().size());
}

const std::vector<BaseNode::WPtr>& BaseNode::adjacent_nodes() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::adjacent_nodes: called for not finalized BaseNode");
    }
    return m_adjacent_nodes;
}

const std::vector<Block::WPtr> &BaseNode::adjacent_blocks() const {
    if (m_editable) {
        throw std::runtime_error("BaseNode::adjacent_blocks: called for not finalized BaseNode");
    }
    return m_adjacent_blocks;
}

} // namespace zephyr::geom::generator
