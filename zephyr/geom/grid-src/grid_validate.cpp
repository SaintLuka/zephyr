#include <zephyr/geom/grid.h>

namespace zephyr::geom {

namespace {

void append_line(std::string* out, const std::string& s) {
    if (!out) return;
    if (!out->empty()) out->push_back('\n');
    *out += s;
}

bool is_finite(double x) {
    return std::isfinite(x);
}

bool is_finite(const Vector3d& v) {
    return is_finite(v.x()) && is_finite(v.y()) && is_finite(v.z());
}

// ------------------------ CellType arity helpers ------------------------
bool cell_nodes_count_ok(CellType t, std::size_t n) {
    switch (t) {
        case CellType::TRIANGLE: return n == 3;
        case CellType::QUAD:     return n == 4;
        case CellType::POLYGON:  return n >= 3;
        default:
            // Сообщим в report, что тип неизвестен (см. ниже).
            return n > 0;
    }
}

// Число гран условий
bool cell_face_bc_count_ok(CellType t, std::size_t node_count, std::size_t bc_count) {
    if (bc_count == 0) return true; // разрешаем "не задано"
    switch (t) {
        case CellType::TRIANGLE: return bc_count == 3;
        case CellType::QUAD:     return bc_count == 4;
        case CellType::POLYGON:  return bc_count == node_count;
        default:
            return true;
    }
}

// ---------- CSR basic checks (предполагаем семантику csr.size() == n_rows()) ----------
template <class CsrT>
bool validate_csr_shape(const CsrT& csr, std::size_t expected_rows, std::string* report, const char* name) {
    const auto& off = csr.offsets;
    const auto& val = csr.values;

    bool ok = true;

    if (csr.size() != expected_rows) {
        ok = false;
        append_line(report, std::string(name) + ": csr.size() != expected_rows");
    }

    if (off.size() != expected_rows + 1) {
        ok = false;
        append_line(report, std::string(name) + ": offsets.size() != rows+1");
        return ok; // дальнейшие проверки бессмысленны
    }

    if (!off.empty() && off.front() != 0) {
        ok = false;
        append_line(report, std::string(name) + ": offsets[0] must be 0");
    }

    for (std::size_t i = 0; i + 1 < off.size(); ++i) {
        if (off[i] > off[i + 1]) {
            ok = false;
            append_line(report, std::string(name) + ": offsets must be non-decreasing");
            break;
        }
    }

    if (!off.empty() && static_cast<std::size_t>(off.back()) != val.size()) {
        ok = false;
        append_line(report, std::string(name) + ": offsets.back() must equal values.size()");
    }

    return ok;
}

template <typename RowView>
bool contains_id(const RowView& row, id_t x) {
    for (const auto& v : row) {
        if (v == x) return true;
    }
    return false;
}

} // anonymous namespace

// ============================================================================
//                 Grid validate: level 0 default validate()
// ===========================================================================

bool Grid::validate(std::string* report) const {
    if (m_state == State::editable) {
        return validate_draft(report);
    }
    return validate_finalized_basic(report);
}

// ============================================================================
//                 Grid validate: Level 1 (draft / pre-finalize)
// ============================================================================

bool Grid::validate_draft(std::string* report) const {
    if (report) report->clear();

    if (m_state != State::editable) {
        append_line(report, "validate_draft: grid is not editable.");
        return false;
    }

    if (!m_draft.has_value()) {
        append_line(report, "validate_draft: m_draft is empty.");
        return false;
    }

    const DraftData& d = *m_draft;

    bool ok = true;

    // unique_nodes pointers + builder state sanity
    for (id_t i = 0; i < d.nodes.size(); ++i) {
        auto node = d.nodes[i];
        if (!node) {
            ok = false;
            append_line(report, "validate_draft: unique_nodes contains null pointer at index " + std::to_string(i) + ".");
            continue;
        }

        if (node->id() != i) {
            ok = false;
            append_line(report, "validate_draft: NodeInput builder state mismatch at unique_nodes[" + std::to_string(i) + "].");
        }

        if (!is_finite(node->pos)) {
            ok = false;
            append_line(report, "validate_draft: NodeInput position is not finite at unique_nodes[" + std::to_string(i) + "].");
        }
    }

    // draft cells
    for (std::size_t ci = 0; ci < d.cells.size(); ++ci) {
        const DraftCell& c = d.cells[ci];

        const std::size_t n = c.node_ids.size();
        if (!cell_nodes_count_ok(c.type, n)) {
            ok = false;
            append_line(report,
                "validate_draft: cell " + std::to_string(ci) +
                " has invalid node count (" + std::to_string(n) + ").");
        }

        // node id bounds
        for (std::size_t k = 0; k < c.node_ids.size(); ++k) {
            const id_t nid = c.node_ids[k];
            if (nid == invalid_id || static_cast<std::size_t>(nid) >= d.nodes.size()) {
                ok = false;
                append_line(report,
                    "validate_draft: cell " + std::to_string(ci) +
                    " has out-of-range node id at k=" + std::to_string(k) + ".");
            }
        }

        // face_bc sizing (for types where you currently pass it)
        if (!cell_face_bc_count_ok(c.type, c.node_ids.size(), c.face_bc.size())) {
            ok = false;
            append_line(report,
                "validate_draft: cell " + std::to_string(ci) +
                " face_bc size mismatch (face_bc=" + std::to_string(c.face_bc.size()) +
                ", nodes=" + std::to_string(c.node_ids.size()) + ").");
        }

        // duplicate nodes inside a cell (debug-ish, но полезно ловит ошибки на входе)
        {
            std::vector<id_t> tmp = c.node_ids;
            std::ranges::sort(tmp);
            if (std::ranges::adjacent_find(tmp) != tmp.end()) {
                ok = false;
                append_line(report, "validate_draft: cell " + std::to_string(ci) + " contains duplicate node ids.");
            }
        }
    }

    return ok;
}

// ============================================================================
//    Grid validate: Level 2 (finalized basic: sizes + options compliance)
// ============================================================================

bool Grid::validate_finalized_basic(std::string* report) const {
    if (report) report->clear();

    if (m_state != State::finalized) {
        append_line(report, "validate_finalized_basic: grid is not finalized.");
        return false;
    }

    const BuildOptions& opt = m_options;

    // validate option compatibility
    try {
        opt.validate_or_throw(m_dim);
    } catch (const std::exception& e) {
        append_line(report, std::string("validate_finalized_basic: invalid BuildOptions: ") + e.what());
        return false;
    }

    bool ok = true;

    if (m_dim < 2) {
        ok = false;
        append_line(report, "validate_finalized_basic: strange dimension " + std::to_string(m_dim) + ".");
    }
    if (m_type == Type::NONE) {
        ok = false;
        append_line(report, "validate_finalized_basic: strange grid type.");
    }

    // Nodes basic
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {
        if (!is_finite(m_nodes[i].pos)) {
            ok = false;
            append_line(report, "validate_finalized_basic: node " + std::to_string(i) + " position is not finite.");
        }
    }

    // Cells basic: bounds + arity
    for (std::size_t ci = 0; ci < m_cells.size(); ++ci) {
        const Cell& c = m_cells[ci];
        const std::size_t n = c.nodes.size();

        if (!cell_nodes_count_ok(c.type, n)) {
            ok = false;
            append_line(report,
                "validate_finalized_basic: cell " + std::to_string(ci) +
                " has invalid node count (" + std::to_string(n) + ").");
        }

        for (std::size_t k = 0; k < c.nodes.size(); ++k) {
            const id_t nid = c.nodes[k];
            if (nid == invalid_id || static_cast<std::size_t>(nid) >= m_nodes.size()) {
                ok = false;
                append_line(report,
                    "validate_finalized_basic: cell " + std::to_string(ci) +
                    " has out-of-range node id at k=" + std::to_string(k) + ".");
            }
        }
    }

    // Faces cache compliance
    if (opt.faces == BuildOptions::FaceOption::none) {
        if (m_faces.has_value()) {
            ok = false;
            append_line(report, "validate_finalized_basic: opt.faces==none but m_faces is present.");
        }
    } else {
        if (!m_faces.has_value()) {
            ok = false;
            append_line(report, "validate_finalized_basic: opt.faces!=none but m_faces is missing.");
        } else {
            // structural shape
            try {
                m_faces->validate_or_throw(m_cells.size());
            } catch (const std::exception& e) {
                ok = false;
                append_line(report, std::string("validate_finalized_basic: FacesCache invalid: ") + e.what());
            }

            // option-specific presence
            if (opt.build_twin_face && m_faces->twin_face.empty()) {
                ok = false;
                append_line(report, "validate_finalized_basic: build_twin_face requested but twin_face is empty.");
            }
            if (opt.build_face_local_indices && m_faces->local_nodes.size() != m_faces->n_faces()) {
                ok = false;
                append_line(report, "validate_finalized_basic: build_face_local_indices requested but local_nodes not built.");
            }
        }
    }

    // Edges cache compliance
    if (!opt.build_edges) {
        if (m_edges.has_value()) {
            ok = false;
            append_line(report, "validate_finalized_basic: opt.build_edges==false but m_edges is present.");
        }
    }
    else {
        if (!m_edges.has_value()) {
            ok = false;
            append_line(report, "validate_finalized_basic: opt.build_edges=true but m_edges is missing.");
        } else {
            // structural shape
            try {
                m_edges->validate_or_throw();
            } catch (const std::exception& e) {
                ok = false;
                append_line(report, std::string("validate_finalized_basic: EdgesCache invalid: ") + e.what());
            }
        }
    }

    // Geometry cache compliance
    if (!m_geom.has_value()) {
        // allow empty geometry always
    } else {
        const GeometryCache& g = *m_geom;
        const std::size_t nc = m_cells.size();

        if (!g.cell_centroids.empty() && g.cell_centroids.size() != nc) {
            ok = false;
            append_line(report, "validate_finalized_basic: cell_centroids size != n_cells.");
        }
        if (!g.cell_volumes.empty() && g.cell_volumes.size() != nc) {
            ok = false;
            append_line(report, "validate_finalized_basic: cell_volumes size != n_cells.");
        }

        const bool face_geom_present = !g.face_areas.empty() || !g.face_centroids.empty() || !g.face_normals.empty();
        if (face_geom_present) {
            if (!m_faces.has_value()) {
                ok = false;
                append_line(report, "validate_finalized_basic: face geometry present but faces cache is missing.");
            } else {
                const std::size_t nf = m_faces->n_faces();
                if (!g.face_areas.empty() && g.face_areas.size() != nf) {
                    ok = false;
                    append_line(report, "validate_finalized_basic: face_areas size != n_faces.");
                }
                if (!g.face_centroids.empty() && g.face_centroids.size() != nf) {
                    ok = false;
                    append_line(report, "validate_finalized_basic: face_centroids size != n_faces.");
                }
                if (!g.face_normals.empty() && g.face_normals.size() != nf) {
                    ok = false;
                    append_line(report, "validate_finalized_basic: face_normals size != n_faces.");
                }
            }
        }

        // if options asked for geometry but cache doesn't contain it
        if (opt.compute_cell_geometry) {
            if (g.cell_centroids.empty() || g.cell_volumes.empty()) {
                ok = false;
                append_line(report, "validate_finalized_basic: compute_cell_geometry requested but cell geometry arrays are empty.");
            }
        }
        if (opt.compute_face_geometry) {
            if (g.face_centroids.empty() || g.face_normals.empty() || g.face_areas.empty()) {
                ok = false;
                append_line(report, "validate_finalized_basic: compute_face_geometry requested but face geometry arrays are empty.");
            }
        }
    }

    // Incidence cache compliance
    if (opt.build_node_cells || opt.build_node_faces) {
        if (!m_inc.has_value()) {
            ok = false;
            append_line(report, "validate_finalized_basic: incidence requested but m_inc is missing.");
        } else {
            const std::size_t nn = m_nodes.size();

            ok = validate_csr_shape(m_inc->node_cells, nn, report, "IncidenceCache.node_cells") && ok;

            if (opt.build_node_faces) {
                if (!m_faces.has_value()) {
                    ok = false;
                    append_line(report, "validate_finalized_basic: build_node_faces requested but faces cache missing.");
                }
                ok = validate_csr_shape(m_inc->node_faces, nn, report, "IncidenceCache.node_faces") && ok;
            }
        }
    } else {
        // allowed: m_inc absent or present-but-empty; choose your policy.
    }

    return ok;
}

// ============================================================================
//           Grid validate: Level 3 (finalized full debug checkout)
// ============================================================================

bool Grid::validate_finalized_full(std::string* report) const {
    if (report) report->clear();

    // basic first
    std::string tmp;
    if (!validate_finalized_basic(&tmp)) {
        append_line(report, "validate_finalized_full: basic validation failed:");
        append_line(report, tmp);
        return false;
    }

    bool ok = true;

    if (m_type == Type::NONE) {
        ok = false;
        append_line(report, "validate_finalized_full: wrong grid type.");
    }

    // Cells, dimension, type
    for (std::size_t ci = 0; ci < m_cells.size(); ++ci) {
        auto ctype = m_cells[ci].type;
        if (get_dimension(ctype) != m_dim) {
            ok = false;
            append_line(report, "validate_finalized_full: cell " + std::to_string(ci) + " wrong dimension.");
        }
        bool bad_type = false;
        if (m_type == Type::AMR) {
            if (ctype != CellType::AMR2D && ctype != CellType::AMR3D)
                bad_type = true;
        }
        else if (m_type == Type::TRI) {
            if (ctype != CellType::TRIANGLE && ctype != CellType::TETRA) {
                bad_type = true;
            }
        }
        else if (m_type == Type::QUAD) {
            if (ctype != CellType::QUAD && ctype != CellType::HEXAHEDRON) {
                bad_type = true;
            }
        }
        else if (m_type == Type::POLY) {
            if (ctype == CellType::AMR2D || ctype == CellType::AMR3D || ctype == CellType::POLYHEDRON) {
                bad_type = true;
            }
        }
        else if (m_type == Type::POLYHEDRON) {
            if (ctype == CellType::AMR2D || ctype == CellType::AMR3D) {
                bad_type = true;
            }
        }
        if (bad_type) {
            ok = false;
            append_line(report, "validate_finalized_full: cell " + std::to_string(ci) + " grid type mismatch.");
        }
    }

    // Cells: duplicates, finiteness, etc.
    for (std::size_t ci = 0; ci < m_cells.size(); ++ci) {
        const Cell& c = m_cells[ci];

        std::vector<id_t> tmp_nodes(c.nodes.begin(), c.nodes.end());
        std::sort(tmp_nodes.begin(), tmp_nodes.end());
        if (std::adjacent_find(tmp_nodes.begin(), tmp_nodes.end()) != tmp_nodes.end()) {
            ok = false;
            append_line(report, "validate_finalized_full: cell " + std::to_string(ci) + " contains duplicate node ids.");
        }
    }

    // Geometry deep sanity
    if (m_geom.has_value()) {
        const GeometryCache& g = *m_geom;

        for (std::size_t i = 0; i < g.cell_centroids.size(); ++i) {
            if (!is_finite(g.cell_centroids[i])) {
                ok = false;
                append_line(report, "validate_finalized_full: cell centroid not finite at cell " + std::to_string(i) + ".");
            }
        }
        for (std::size_t i = 0; i < g.cell_volumes.size(); ++i) {
            if (!is_finite(g.cell_volumes[i])) {
                ok = false;
                append_line(report, "validate_finalized_full: cell volume not finite at cell " + std::to_string(i) + ".");
            }
        }

        for (std::size_t i = 0; i < g.face_centroids.size(); ++i) {
            if (!is_finite(g.face_centroids[i])) {
                ok = false;
                append_line(report, "validate_finalized_full: face centroid not finite at face " + std::to_string(i) + ".");
            }
        }
        for (std::size_t i = 0; i < g.face_normals.size(); ++i) {
            if (!is_finite(g.face_normals[i])) {
                ok = false;
                append_line(report, "validate_finalized_full: face normal not finite at face " + std::to_string(i) + ".");
            }
        }
        for (std::size_t i = 0; i < g.face_areas.size(); ++i) {
            if (!is_finite(g.face_areas[i])) {
                ok = false;
                append_line(report, "validate_finalized_full: face area not finite at face " + std::to_string(i) + ".");
            }
        }
    }

    const BuildOptions& opt = m_options;

    // Faces deep checks (если faces строились)
    if (opt.faces != BuildOptions::FaceOption::none && m_faces.has_value()) {
        const FacesCache& f = *m_faces;
        const std::size_t nf = f.n_faces();
        const std::size_t nc = m_cells.size();
        const std::size_t nn = m_nodes.size();

        // Owner/neighbor bounds + face_nodes node id bounds
        for (std::size_t face_id = 0; face_id < nf; ++face_id) {
            const id_t owner = f.owner_cell[face_id];
            const id_t neigh = f.neighbor_cell[face_id];

            if (owner == invalid_id || static_cast<std::size_t>(owner) >= nc) {
                ok = false;
                append_line(report, "validate_finalized_full: face " + std::to_string(face_id) + " owner out of range.");
                continue; // дальше owner нужен
            }
            if (neigh != invalid_id && static_cast<std::size_t>(neigh) >= nc) {
                ok = false;
                append_line(report, "validate_finalized_full: face " + std::to_string(face_id) + " neighbor out of range.");
            }
            if (neigh == owner) {
                ok = false;
                append_line(report, "validate_finalized_full: face " + std::to_string(face_id) + " neighbor == owner.");
            }

            // face must appear in owner cell_faces row
            if (!contains_id(f.cell_faces[owner], static_cast<id_t>(face_id))) {
                ok = false;
                append_line(report, "validate_finalized_full: face " + std::to_string(face_id) + " not found in cell_faces[owner].");
            }

            // In unique mode: face should appear in neighbor row too (for internal face)
            if (opt.faces == BuildOptions::FaceOption::unique && neigh != invalid_id) {
                if (!contains_id(f.cell_faces[neigh], static_cast<id_t>(face_id))) {
                    ok = false;
                    append_line(report, "validate_finalized_full: (unique) face " + std::to_string(face_id) + " not found in cell_faces[neighbor].");
                }
            }

            // face_nodes validity
            for (id_t nid : f.face_nodes[face_id]) {
                if (nid == invalid_id || static_cast<std::size_t>(nid) >= nn) {
                    ok = false;
                    append_line(report, "validate_finalized_full: face " + std::to_string(face_id) + " has out-of-range node id in face_nodes.");
                    break;
                }
            }
        }

        // twin_face symmetry (если построен)
        if (opt.build_twin_face && !f.twin_face.empty()) {
            for (std::size_t face_id = 0; face_id < nf; ++face_id) {
                const id_t twin = f.twin_face[face_id];
                if (twin == invalid_id) continue;
                if (static_cast<std::size_t>(twin) >= nf) {
                    ok = false;
                    append_line(report, "validate_finalized_full: twin_face out of range at face " + std::to_string(face_id) + ".");
                    continue;
                }
                if (f.twin_face[twin] != static_cast<id_t>(face_id)) {
                    ok = false;
                    append_line(report, "validate_finalized_full: twin symmetry broken at face " + std::to_string(face_id) + ".");
                }
                // owner/neigh should be swapped (common invariant for per_cell)
                if (!(f.owner_cell[face_id] == f.neighbor_cell[twin] && f.neighbor_cell[face_id] == f.owner_cell[twin])) {
                    ok = false;
                    append_line(report, "validate_finalized_full: twin owner/neighbor swap invariant broken at face " + std::to_string(face_id) + ".");
                }
            }
        }

        // local_nodes correctness (если построен)
        if (opt.build_face_local_indices && f.local_nodes.size() == nf) {
            for (std::size_t face_id = 0; face_id < nf; ++face_id) {
                const id_t owner = f.owner_cell[face_id];
                if (owner == invalid_id || static_cast<std::size_t>(owner) >= nc) continue;

                const auto face_nodes_row = f.face_nodes[face_id];
                const auto local_row      = f.local_nodes[face_id];

                // должны быть одинаковой длины
                if (face_nodes_row.size() != local_row.size()) {
                    ok = false;
                    append_line(report, "validate_finalized_full: local_nodes row size != face_nodes row size at face " + std::to_string(face_id) + ".");
                    continue;
                }

                const Cell& cell = m_cells[owner];

                std::size_t k = 0;
                auto itN = face_nodes_row.begin();
                auto itL = local_row.begin();
                for (; itN != face_nodes_row.end() && itL != local_row.end(); ++itN, ++itL, ++k) {
                    const id_t global_n = *itN;
                    const std::uint16_t local_i = *itL;

                    if (static_cast<std::size_t>(local_i) >= cell.nodes.size()) {
                        ok = false;
                        append_line(report,
                            "validate_finalized_full: local index out of range at face " + std::to_string(face_id) +
                            ", k=" + std::to_string(k) + ".");
                        break;
                    }

                    if (cell.nodes[local_i] != global_n) {
                        ok = false;
                        append_line(report,
                            "validate_finalized_full: local_nodes mismatch at face " + std::to_string(face_id) +
                            ", k=" + std::to_string(k) + ".");
                        break;
                    }
                }
            }
        }
    }

    if (m_edges.has_value()) {
        // edges meaningful only for 3D
        if (m_dim != 3) {
            ok = false;
            append_line(report, "validate_finalized_full: EdgesCache present but grid dimension != 3.");
        }

        const EdgesCache& e = *m_edges;
        const std::size_t ne = e.n_edges();
        const std::size_t nn = m_nodes.size();
        const std::size_t nc = m_cells.size();

        // Базовая форма/согласованность самого кэша
        try {
            e.validate_or_throw();
        } catch (const std::exception& ex) {
            ok = false;
            append_line(report, std::string("validate_finalized_full: EdgesCache invalid: ") + ex.what());
        }

        // Проверка индексов узлов ребра
        for (std::size_t ei = 0; ei < ne; ++ei) {
            const id_t a = e.edges[ei][0];
            const id_t b = e.edges[ei][1];

            if (a == invalid_id || b == invalid_id ||
                static_cast<std::size_t>(a) >= nn || static_cast<std::size_t>(b) >= nn) {
                ok = false;
                append_line(report, "validate_finalized_full: edge node id out of range at edge " + std::to_string(ei) + ".");
                continue;
                }
            if (a == b) {
                ok = false;
                append_line(report, "validate_finalized_full: degenerate edge (a==b) at edge " + std::to_string(ei) + ".");
            }
        }

        // Проверка CSR shape (если заполнены)
        // edge_faces: rows == n_edges, entries are face_ids
        if (!e.edge_faces.offsets.empty() || !e.edge_faces.values.empty()) {
            ok = validate_csr_shape(e.edge_faces, ne, report, "EdgesCache.edge_faces") && ok;

            // если есть faces cache — проверяем диапазоны face_id
            if (m_faces.has_value()) {
                const std::size_t nf = m_faces->n_faces();
                for (std::size_t ei = 0; ei < ne; ++ei) {
                    for (id_t fid : e.edge_faces[ei]) {
                        if (fid == invalid_id || static_cast<std::size_t>(fid) >= nf) {
                            ok = false;
                            append_line(report, "validate_finalized_full: edge_faces contains out-of-range face_id at edge "
                                                + std::to_string(ei) + ".");
                            break;
                        }
                    }
                }
            }
        }

        // edge_cells: rows == n_edges, entries are cell_ids
        if (!e.edge_cells.offsets.empty() || !e.edge_cells.values.empty()) {
            ok = validate_csr_shape(e.edge_cells, ne, report, "EdgesCache.edge_cells") && ok;

            for (std::size_t ei = 0; ei < ne; ++ei) {
                for (id_t cid : e.edge_cells[ei]) {
                    if (cid == invalid_id || static_cast<std::size_t>(cid) >= nc) {
                        ok = false;
                        append_line(report, "validate_finalized_full: edge_cells contains out-of-range cell_id at edge "
                                            + std::to_string(ei) + ".");
                        break;
                    }
                }
            }
        }

        // 3) Согласованность с ячейками/гранями (дорого, но это full-check)
        // Включаем только если есть faces cache (для проверки edge_faces) и/или cells (для edge_cells).
        auto cell_contains_nodes = [&](id_t cid, id_t a, id_t b) -> bool {
            const Cell& c = m_cells[static_cast<std::size_t>(cid)];
            bool ha = false, hb = false;
            for (id_t nid : c.nodes) {
                ha = ha || (nid == a);
                hb = hb || (nid == b);
                if (ha && hb) return true;
            }
            return false;
        };

        auto face_contains_nodes = [&](id_t fid, id_t a, id_t b) -> bool {
            // require faces cache
            const auto row = m_faces->face_nodes[static_cast<std::size_t>(fid)];
            bool ha = false, hb = false;
            for (id_t nid : row) {
                ha = ha || (nid == a);
                hb = hb || (nid == b);
                if (ha && hb) return true;
            }
            return false;
        };

        // edge_cells: каждый cell в строке должен реально содержать оба узла ребра
        if (!e.edge_cells.offsets.empty()) {
            for (std::size_t ei = 0; ei < ne; ++ei) {
                const id_t a = e.edges[ei][0];
                const id_t b = e.edges[ei][1];
                for (id_t cid : e.edge_cells[ei]) {
                    if (cid == invalid_id || static_cast<std::size_t>(cid) >= nc) continue;
                    if (!cell_contains_nodes(cid, a, b)) {
                        ok = false;
                        append_line(report, "validate_finalized_full: edge_cells inconsistent: cell does not contain edge nodes "
                                            "(edge " + std::to_string(ei) + ", cell " + std::to_string(cid) + ").");
                        break;
                    }
                }
            }
        }

        // edge_faces: каждый face в строке должен реально содержать оба узла ребра
        if (!e.edge_faces.offsets.empty()) {
            if (!m_faces.has_value()) {
                ok = false;
                append_line(report, "validate_finalized_full: edge_faces present but FacesCache is missing.");
            } else {
                const std::size_t nf = m_faces->n_faces();
                for (std::size_t ei = 0; ei < ne; ++ei) {
                    const id_t a = e.edges[ei][0];
                    const id_t b = e.edges[ei][1];
                    for (id_t fid : e.edge_faces[ei]) {
                        if (fid == invalid_id || static_cast<std::size_t>(fid) >= nf) continue;
                        if (!face_contains_nodes(fid, a, b)) {
                            ok = false;
                            append_line(report, "validate_finalized_full: edge_faces inconsistent: face does not contain edge nodes "
                                                "(edge " + std::to_string(ei) + ", face " + std::to_string(fid) + ").");
                            break;
                        }
                    }
                }
            }
        }
    }

    // Incidence deep checks (если построены)
    if (m_inc.has_value()) {
        const IncidenceCache& inc = *m_inc;

        // node_cells: каждый cell должен присутствовать в строках своих узлов
        if (opt.build_node_cells) {
            for (std::size_t ci = 0; ci < m_cells.size(); ++ci) {
                const Cell& c = m_cells[ci];
                for (id_t nid : c.nodes) {
                    if (!contains_id(inc.node_cells[nid], static_cast<id_t>(ci))) {
                        ok = false;
                        append_line(report, "validate_finalized_full: node_cells missing (node " + std::to_string(nid) +
                                             " does not contain cell " + std::to_string(ci) + ").");
                    }
                }
            }
        }

        // node_faces: если есть faces
        if (opt.build_node_faces && m_faces.has_value()) {
            const FacesCache& f = *m_faces;
            for (std::size_t face_id = 0; face_id < f.n_faces(); ++face_id) {
                for (id_t nid : f.face_nodes[face_id]) {
                    if (!contains_id(inc.node_faces[nid], static_cast<id_t>(face_id))) {
                        ok = false;
                        append_line(report, "validate_finalized_full: node_faces missing (node " + std::to_string(nid) +
                                             " does not contain face " + std::to_string(face_id) + ").");
                    }
                }
            }
        }
    }

    return ok;
}

} // namespace zephyr::geom