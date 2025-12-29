#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cstdint>

#include <zephyr/geom/grid.h>

namespace zephyr::geom {

void Grid::finalize(const BuildOptions& opt) {
    require_editable_();
    opt.validate_or_throw(m_dim);
    m_options = opt; // save build options

    if (!m_draft) {
        throw std::runtime_error("Grid::finalize: no draft data");
    }

    // Validate draft
    std::string report;
    if (bool valid = validate_draft(&report); !valid) {
        std::cerr << report << std::endl;
        throw std::runtime_error("Grid::finalize(): Bad draft");
    }

    DraftData& d = *m_draft;
    const id_t n_cells = static_cast<id_t>(d.cells.size());
    const id_t n_nodes = static_cast<id_t>(d.nodes.size());

    // -----------------------------
    // helpers
    // -----------------------------
    auto throw_if = [](bool cond, const char* msg) {
        if (cond) throw std::runtime_error(msg);
    };

    auto node_pos = [&](id_t nid) -> const Vector3d& {
        return m_nodes[nid].pos;
    };

    auto mean_point = [&](const std::vector<id_t>& ids) -> Vector3d {
        Vector3d c(0,0,0);
        for (id_t id : ids) c += node_pos(id);
        if (!ids.empty()) c /= double(ids.size());
        return c;
    };

    auto safe_normalize = [&](const Vector3d& v) -> Vector3d {
        const double n = v.norm();
        if (n <= 0.0) return Vector3d(0,0,0);
        return v / n;
    };

    // corners for common high-order faces: 3x3 quad face
    auto corners_from_face_nodes = [&](const std::vector<id_t>& fnodes) -> std::vector<id_t> {
        if (fnodes.size() == 9) {
            // row-major 3x3: (0,0)=0 (2,0)=2 (2,2)=8 (0,2)=6
            return { fnodes[0], fnodes[2], fnodes[8], fnodes[6] };
        }
        return fnodes; // assume boundary order already
    };

    // flip orientation of a face node list.
    // For 9-node (3x3) face stored row-major: flip by reversing columns in each row.
    auto flip_face_orientation_inplace = [&](std::vector<id_t>& fnodes) {
        if (fnodes.size() == 9) {
            // rows: 0..2, cols: 0..2
            for (int r = 0; r < 3; ++r) {
                std::swap(fnodes[r*3 + 0], fnodes[r*3 + 2]);
            }
        } else {
            std::reverse(fnodes.begin(), fnodes.end());
        }
    };

    // compute face normal+area+centroid from CORNER cycle (or whole list if it is boundary)
    auto face_geom_from_nodes = [&](const std::vector<id_t>& fnodes_cycle,
                                   const Vector3d& cell_centroid_guess) -> std::tuple<Vector3d,double,Vector3d>
    {
        const std::vector<id_t> cyc = corners_from_face_nodes(fnodes_cycle);
        throw_if(cyc.size() < 2, "Grid::finalize: face has <2 nodes");

        // centroid (simple average of the corner nodes, ok for planar faces)
        Vector3d fc = mean_point(cyc);

        if (m_dim == 2) {
            // 2D "face" is an edge; treat "area" as length, normal in XY plane.
            const Vector3d p0 = node_pos(cyc.front());
            const Vector3d p1 = node_pos(cyc.back());
            const Vector3d t  = (p1 - p0);
            const double len  = t.norm();

            Vector3d n = Vector3d(0,0,0);
            if (len > 0.0) {
                const Vector3d z(0,0,1);
                // For CCW boundary order, outward ~ t x z (see note in analysis).
                n = safe_normalize(t.cross(z));
                if (n.dot(fc - cell_centroid_guess) < 0.0) n = -n;
            }
            return { n, len, fc };
        }

        // m_dim == 3: Newell normal on polygon corners
        Vector3d nsum(0,0,0);
        for (size_t i = 0; i < cyc.size(); ++i) {
            const Vector3d p  = node_pos(cyc[i]);
            const Vector3d q  = node_pos(cyc[(i+1) % cyc.size()]);
            nsum.x() += (p.y() - q.y()) * (p.z() + q.z());
            nsum.y() += (p.z() - q.z()) * (p.x() + q.x());
            nsum.z() += (p.x() - q.x()) * (p.y() + q.y());
        }
        const double nn = nsum.norm();
        double area = 0.5 * nn;
        Vector3d n = (nn > 0.0) ? (nsum / nn) : Vector3d(0,0,0);

        // ensure outward vs cell_centroid_guess
        if (n.dot(fc - cell_centroid_guess) < 0.0) n = -n;

        return { n, area, fc };
    };

    // compute cell centroid guess from its nodes (cheap)
    auto cell_centroid_guess = [&](id_t cid) -> Vector3d {
        const Cell& c = m_cells[cid];
        Vector3d cc(0,0,0);
        for (id_t nid : c.nodes) cc += node_pos(nid);
        if (!c.nodes.empty()) cc /= double(c.nodes.size());
        return cc;
    };

    // CSR builder from vector<vector<id_t>>
    auto build_csr = [&](Csr<id_t, id_t>& csr, const std::vector<std::vector<id_t>>& rows) {
        csr.offsets.assign(rows.size() + 1, 0);
        for (size_t r = 0; r < rows.size(); ++r) {
            csr.offsets[r+1] = csr.offsets[r] + static_cast<id_t>(rows[r].size());
        }
        csr.values.resize(csr.offsets.back());
        id_t k = 0;
        for (size_t r = 0; r < rows.size(); ++r) {
            for (id_t v : rows[r]) csr.values[k++] = v;
        }
    };

    // Node->(cells/faces) CSR builder (two-pass)
    auto build_node_incidence = [&](Csr<id_t,id_t>& out,
                                   id_t rows_count,
                                   auto foreach_relation /* fn(push(row, value)) */)
    {
        std::vector<id_t> counts(rows_count, 0);
        foreach_relation([&](id_t row, id_t /*val*/){ ++counts[row]; });

        out.offsets.assign(rows_count + 1, 0);
        for (id_t i = 0; i < rows_count; ++i) out.offsets[i+1] = out.offsets[i] + counts[i];

        out.values.assign(out.offsets.back(), invalid_id);
        std::vector<id_t> cursor = out.offsets;

        foreach_relation([&](id_t row, id_t val){
            out.values[cursor[row]++] = val;
        });
    };

    // FaceKey for merging (order-insensitive: sorted node ids)
    struct FaceKey {
        std::vector<id_t> ids;
        bool operator==(const FaceKey& o) const noexcept { return ids == o.ids; }
    };
    struct FaceKeyHash {
        std::size_t operator()(const FaceKey& k) const noexcept {
            // FNV-1a-ish combine
            std::size_t h = 1469598103934665603ull;
            for (id_t v : k.ids) {
                h ^= std::size_t(v) + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
            }
            return h;
        }
    };
    auto make_face_key = [&](const std::vector<id_t>& fnodes) -> FaceKey {
        FaceKey k{fnodes};
        std::sort(k.ids.begin(), k.ids.end());
        k.ids.erase(std::unique(k.ids.begin(), k.ids.end()), k.ids.end());
        return k;
    };

    // local indices in cell.nodes (uint8)
    auto build_local_indices = [&](id_t cid, const std::vector<id_t>& fnodes) -> std::vector<std::uint8_t> {
        std::vector<std::uint8_t> loc;
        loc.reserve(fnodes.size());
        const Cell& c = m_cells[cid];
        for (id_t nid : fnodes) {
            int found = -1;
            for (int i = 0; i < (int)c.nodes.size(); ++i) {
                if (c.nodes[i] == nid) { found = i; break; }
            }
            throw_if(found < 0 || found > 255, "Grid::finalize: local face node index overflow/invalid");
            loc.push_back(static_cast<std::uint8_t>(found));
        }
        return loc;
    };

    // -----------------------------
    // 1) build base nodes/cells
    // -----------------------------
    m_nodes.clear();
    m_cells.clear();
    m_nodes.resize(n_nodes);
    m_cells.resize(n_cells);

    for (id_t i = 0; i < n_nodes; ++i) {
        const NodeInput* in = d.nodes[i];
        m_nodes[i].pos = in->pos;
        m_nodes[i].bc  = in->bc;
    }

    for (id_t cid = 0; cid < n_cells; ++cid) {
        const DraftCell& dc = d.cells[cid];
        Cell& c = m_cells[cid];
        c.type = dc.type;
        c.nodes.clear();
        c.nodes.insert(c.nodes.end(), dc.node_ids.begin(), dc.node_ids.end());
        throw_if(c.nodes.empty(), "Grid::finalize: cell has no nodes");
    }

    // clear old caches
    m_faces.reset();
    m_edges.reset();
    m_inc.reset();
    m_geom.reset();

    // -----------------------------
    // 2) prepare per-cell temporary faces (needed for:
    //    - building FacesCache
    //    - computing cell geometry robustly even in unique mode
    // -----------------------------
    struct TempFace {
        std::vector<id_t> nodes;                  // oriented for THIS cell
        std::vector<std::uint8_t> local_indices;  // optional
        Boundary bc{Boundary::UNDEFINED};
    };
    std::vector<std::vector<TempFace>> tmp_cell_faces;
    const bool need_tmp_faces =
        (opt.faces != BuildOptions::FaceOption::none) ||
        opt.compute_cell_geometry ||
        (opt.build_edges); // edges from faces is easiest when possible

    if (need_tmp_faces) {
        tmp_cell_faces.resize(n_cells);

        for (id_t cid = 0; cid < n_cells; ++cid) {
            const DraftCell& dc = d.cells[cid];
            const Cell& c = m_cells[cid];

            auto get_bc = [&](size_t lf) -> Boundary {
                if (dc.face_bc.empty()) return Boundary::UNDEFINED;
                if (lf < dc.face_bc.size()) return dc.face_bc[lf];
                return Boundary::UNDEFINED;
            };

            auto push_face = [&](std::vector<id_t> fn, Boundary bc) {
                // orient via geometric check vs cell centroid guess
                const Vector3d cc = cell_centroid_guess(cid);
                auto [n, a, fc] = face_geom_from_nodes(fn, cc);
                // If degenerate (area/len zero), keep as is.
                if ((m_dim == 3 && a > 0.0) || (m_dim == 2 && a > 0.0)) {
                    if (n.dot(fc - cc) < 0.0) {
                        // In case face_geom already flipped, we keep nodes,
                        // but for consistency force outward by flipping node order:
                        flip_face_orientation_inplace(fn);
                    }
                }

                TempFace tf;
                tf.nodes = std::move(fn);
                tf.bc = bc;
                if (opt.build_face_local_indices) {
                    tf.local_indices = build_local_indices(cid, tf.nodes);
                }
                tmp_cell_faces[cid].push_back(std::move(tf));
            };

            // ---- face extraction by type ----
            switch (dc.type) {
                // 2D
                case CellType::TRIANGLE: {
                    throw_if(c.nodes.size() < 3, "Grid::finalize: TRIANGLE needs 3 nodes");
                    push_face({c.nodes[0], c.nodes[1]}, get_bc(0));
                    push_face({c.nodes[1], c.nodes[2]}, get_bc(1));
                    push_face({c.nodes[2], c.nodes[0]}, get_bc(2));
                } break;

                case CellType::QUAD: {
                    throw_if(c.nodes.size() < 4, "Grid::finalize: QUAD needs 4 nodes");
                    push_face({c.nodes[0], c.nodes[1]}, get_bc(0));
                    push_face({c.nodes[1], c.nodes[2]}, get_bc(1));
                    push_face({c.nodes[2], c.nodes[3]}, get_bc(2));
                    push_face({c.nodes[3], c.nodes[0]}, get_bc(3));
                } break;

                case CellType::POLYGON: {
                    const size_t nn = c.nodes.size();
                    throw_if(nn < 3, "Grid::finalize: POLYGON needs >=3 nodes");
                    for (size_t i = 0; i < nn; ++i) {
                        const id_t a = c.nodes[i];
                        const id_t b = c.nodes[(i+1) % nn];
                        push_face({a, b}, get_bc(i));
                    }
                } break;

                case CellType::AMR2D: {
                    // 9 nodes: 3x3 row-major
                    throw_if(c.nodes.size() != 9, "Grid::finalize: AMR2D expects exactly 9 nodes");
                    // bottom: [0,1,2]
                    push_face({c.nodes[0], c.nodes[1], c.nodes[2]}, get_bc(0));
                    // right: [2,5,8]
                    push_face({c.nodes[2], c.nodes[5], c.nodes[8]}, get_bc(1));
                    // top: [8,7,6]
                    push_face({c.nodes[8], c.nodes[7], c.nodes[6]}, get_bc(2));
                    // left: [6,3,0]
                    push_face({c.nodes[6], c.nodes[3], c.nodes[0]}, get_bc(3));
                } break;

                    // 3D classic
                case CellType::TETRA: {
                    throw_if(c.nodes.size() < 4, "Grid::finalize: TETRA needs 4 nodes");
                    // faces as triangles; orientation fixed later by geometric check
                    push_face({c.nodes[0], c.nodes[2], c.nodes[1]}, get_bc(0));
                    push_face({c.nodes[0], c.nodes[1], c.nodes[3]}, get_bc(1));
                    push_face({c.nodes[1], c.nodes[2], c.nodes[3]}, get_bc(2));
                    push_face({c.nodes[2], c.nodes[0], c.nodes[3]}, get_bc(3));
                } break;

                case CellType::PYRAMID: {
                    throw_if(c.nodes.size() < 5, "Grid::finalize: PYRAMID needs 5 nodes");
                    // base quad + 4 side triangles
                    push_face({c.nodes[0], c.nodes[1], c.nodes[2], c.nodes[3]}, get_bc(0));
                    push_face({c.nodes[0], c.nodes[1], c.nodes[4]}, get_bc(1));
                    push_face({c.nodes[1], c.nodes[2], c.nodes[4]}, get_bc(2));
                    push_face({c.nodes[2], c.nodes[3], c.nodes[4]}, get_bc(3));
                    push_face({c.nodes[3], c.nodes[0], c.nodes[4]}, get_bc(4));
                } break;

                case CellType::WEDGE: {
                    throw_if(c.nodes.size() < 6, "Grid::finalize: WEDGE needs 6 nodes");
                    // 2 triangles + 3 quads
                    push_face({c.nodes[0], c.nodes[2], c.nodes[1]}, get_bc(0));
                    push_face({c.nodes[3], c.nodes[4], c.nodes[5]}, get_bc(1));
                    push_face({c.nodes[0], c.nodes[1], c.nodes[4], c.nodes[3]}, get_bc(2));
                    push_face({c.nodes[1], c.nodes[2], c.nodes[5], c.nodes[4]}, get_bc(3));
                    push_face({c.nodes[2], c.nodes[0], c.nodes[3], c.nodes[5]}, get_bc(4));
                } break;

                case CellType::HEXAHEDRON: {
                    throw_if(c.nodes.size() < 8, "Grid::finalize: HEXAHEDRON needs 8 nodes");
                    // 6 quad faces (orientation fixed by geom check)
                    push_face({c.nodes[0], c.nodes[1], c.nodes[2], c.nodes[3]}, get_bc(0));
                    push_face({c.nodes[4], c.nodes[7], c.nodes[6], c.nodes[5]}, get_bc(1));
                    push_face({c.nodes[0], c.nodes[4], c.nodes[5], c.nodes[1]}, get_bc(2));
                    push_face({c.nodes[1], c.nodes[5], c.nodes[6], c.nodes[2]}, get_bc(3));
                    push_face({c.nodes[2], c.nodes[6], c.nodes[7], c.nodes[3]}, get_bc(4));
                    push_face({c.nodes[3], c.nodes[7], c.nodes[4], c.nodes[0]}, get_bc(5));
                } break;

                case CellType::AMR3D: {
                    throw_if(c.nodes.size() != 27, "Grid::finalize: AMR3D expects exactly 27 nodes");

                    auto idx = [&](int i, int j, int k) -> id_t { return c.nodes[i + 3*j + 9*k]; };

                    auto face9 = [&](auto getter, Boundary bc) {
                        std::vector<id_t> fn;
                        fn.reserve(9);
                        for (int v = 0; v < 3; ++v)
                            for (int u = 0; u < 3; ++u)
                                fn.push_back(getter(u, v));
                        push_face(std::move(fn), bc);
                    };

                    // We generate each face as a 3x3 table (row-major),
                    // then push_face() will flip if needed based on geometry.

                    // z- (k=0): use (u=i, v=j) but with j descending to prefer outward -z
                    face9([&](int u, int v)->id_t {
                        const int i = u;
                        const int j = 2 - v;
                        return idx(i,j,0);
                    }, get_bc(0));

                    // z+ (k=2): (u=i, v=j) j ascending
                    face9([&](int u, int v)->id_t {
                        const int i = u;
                        const int j = v;
                        return idx(i,j,2);
                    }, get_bc(1));

                    // x- (i=0): (u=j, v=k)
                    face9([&](int u, int v)->id_t {
                        const int j = u;
                        const int k = v;
                        return idx(0,j,k);
                    }, get_bc(2));

                    // x+ (i=2): (u=j, v=k)
                    face9([&](int u, int v)->id_t {
                        const int j = u;
                        const int k = v;
                        return idx(2,j,k);
                    }, get_bc(3));

                    // y- (j=0): (u=i, v=k)
                    face9([&](int u, int v)->id_t {
                        const int i = u;
                        const int k = v;
                        return idx(i,0,k);
                    }, get_bc(4));

                    // y+ (j=2): (u=i, v=k)
                    face9([&](int u, int v)->id_t {
                        const int i = u;
                        const int k = v;
                        return idx(i,2,k);
                    }, get_bc(5));

                } break;

                case CellType::POLYHEDRON: {
                    // Draft contains per-face node lists
                    const auto& off = dc.poly_face_offsets;
                    const auto& fn  = dc.poly_face_nodes;
                    throw_if(off.size() < 2, "Grid::finalize: POLYHEDRON offsets empty");
                    const id_t nf = static_cast<id_t>(off.size() - 1);
                    for (id_t f = 0; f < nf; ++f) {
                        const id_t a = off[f];
                        const id_t b = off[f+1];
                        throw_if(b < a, "Grid::finalize: POLYHEDRON bad offsets");
                        std::vector<id_t> face;
                        face.reserve(b - a);
                        for (id_t k = a; k < b; ++k) face.push_back(fn[k]);
                        push_face(std::move(face), get_bc(f));
                    }
                } break;

                default:
                    throw std::runtime_error("Grid::finalize: unsupported cell type in face extraction");
            }
        }
    }

    // -----------------------------
    // 3) build FacesCache (optional)
    // -----------------------------
    if (opt.faces != BuildOptions::FaceOption::none) {
        m_faces.emplace();
        FacesCache& fc = *m_faces;

        std::vector<std::vector<id_t>> cell_face_rows(n_cells);
        std::vector<std::vector<id_t>> face_nodes_rows;
        std::vector<std::vector<std::uint8_t>> face_local_rows;

        fc.owner_cell.clear();
        fc.neighbor_cell.clear();
        fc.face_bc.clear();
        fc.twin_face.clear();
        fc.local_nodes.clear();
        fc.cell_faces.clear();
        fc.face_nodes.clear();

        if (opt.faces == BuildOptions::FaceOption::per_cell) {
            // one face record per cell-face
            const id_t total_faces = [&](){
                id_t s = 0;
                for (id_t cid = 0; cid < n_cells; ++cid) s += (id_t)tmp_cell_faces[cid].size();
                return s;
            }();

            face_nodes_rows.reserve(total_faces);
            if (opt.build_face_local_indices) face_local_rows.reserve(total_faces);

            fc.owner_cell.assign(total_faces, invalid_id);
            fc.neighbor_cell.assign(total_faces, invalid_id);
            fc.face_bc.assign(total_faces, Boundary::UNDEFINED);

            if (opt.build_twin_face) fc.twin_face.assign(total_faces, invalid_id);

            id_t fid = 0;
            for (id_t cid = 0; cid < n_cells; ++cid) {
                cell_face_rows[cid].reserve(tmp_cell_faces[cid].size());
                for (const TempFace& tf : tmp_cell_faces[cid]) {
                    fc.owner_cell[fid] = cid;
                    fc.face_bc[fid]    = tf.bc;

                    face_nodes_rows.push_back(tf.nodes);
                    if (opt.build_face_local_indices) face_local_rows.push_back(tf.local_indices);

                    cell_face_rows[cid].push_back(fid);
                    ++fid;
                }
            }

            // Build adjacency ALWAYS (neighbor_cell always, twin_face optionally)
            {
                std::unordered_map<FaceKey, id_t, FaceKeyHash> map;
                map.reserve((size_t)total_faces * 2);

                for (id_t f = 0; f < total_faces; ++f) {
                    FaceKey key = make_face_key(face_nodes_rows[f]);
                    auto it = map.find(key);
                    if (it == map.end()) {
                        map.emplace(std::move(key), f);
                        continue;
                    }

                    const id_t g = it->second;

                    // Optional: detect non-manifold: more than 2 occurrences
                    // If you want strictness:
                    // if (fc.neighbor_cell[g] != invalid_id) throw std::runtime_error("non-manifold face");

                    // Set neighbors symmetrically
                    fc.neighbor_cell[g] = fc.owner_cell[f];
                    fc.neighbor_cell[f] = fc.owner_cell[g];

                    // Interior => bc undefined (for both records)
                    fc.face_bc[f] = Boundary::UNDEFINED;
                    fc.face_bc[g] = Boundary::UNDEFINED;

                    if (opt.build_twin_face) {
                        fc.twin_face[f] = g;
                        fc.twin_face[g] = f;
                    }
                }
            }
        } else { // unique
            std::unordered_map<FaceKey, id_t, FaceKeyHash> map;
            map.reserve((size_t)n_cells * 8);

            for (id_t cid = 0; cid < n_cells; ++cid) {
                cell_face_rows[cid].reserve(tmp_cell_faces[cid].size());
                for (const TempFace& tf : tmp_cell_faces[cid]) {
                    FaceKey key = make_face_key(tf.nodes);
                    auto it = map.find(key);
                    if (it == map.end()) {
                        const id_t new_id = (id_t)face_nodes_rows.size();
                        map.emplace(std::move(key), new_id);

                        face_nodes_rows.push_back(tf.nodes);
                        if (opt.build_face_local_indices) face_local_rows.push_back(tf.local_indices);

                        fc.owner_cell.push_back(cid);
                        fc.neighbor_cell.push_back(invalid_id);
                        fc.face_bc.push_back(tf.bc);

                        cell_face_rows[cid].push_back(new_id);
                    } else {
                        const id_t fid = it->second;

                        // mark neighbor (non-manifold check)
                        throw_if(fc.owner_cell[fid] == cid, "Grid::finalize: duplicated face in same cell");
                        throw_if(fc.neighbor_cell[fid] != invalid_id && fc.neighbor_cell[fid] != cid,
                                 "Grid::finalize: non-manifold face (more than 2 adjacent cells)");

                        fc.neighbor_cell[fid] = cid;
                        fc.face_bc[fid] = Boundary::UNDEFINED; // interior
                        cell_face_rows[cid].push_back(fid);
                    }
                }
            }
        }

        // pack CSR
        build_csr(fc.cell_faces, cell_face_rows);
        build_csr(fc.face_nodes, face_nodes_rows);

        if (opt.build_face_local_indices) {
            // pack local_nodes CSR<uint8_t>
            fc.local_nodes.offsets.assign(face_local_rows.size() + 1, 0);
            for (size_t r = 0; r < face_local_rows.size(); ++r)
                fc.local_nodes.offsets[r+1] = fc.local_nodes.offsets[r] + (id_t)face_local_rows[r].size();

            fc.local_nodes.values.resize(fc.local_nodes.offsets.back());
            id_t k = 0;
            for (size_t r = 0; r < face_local_rows.size(); ++r) {
                for (auto v : face_local_rows[r]) fc.local_nodes.values[k++] = v;
            }
        }

        fc.validate_or_throw((size_t)n_cells);
    }

    // -----------------------------
    // 4) build IncidenceCache (optional)
    // -----------------------------
    if (opt.build_node_cells || opt.build_node_faces) {
        m_inc.emplace();
        IncidenceCache& ic = *m_inc;
        ic.clear();

        if (opt.build_node_cells) {
            build_node_incidence(
                ic.node_cells, n_nodes,
                [&](auto push){
                    for (id_t cid = 0; cid < n_cells; ++cid) {
                        for (id_t nid : m_cells[cid].nodes) push(nid, cid);
                    }
                }
            );
        }

        if (opt.build_node_faces) {
            throw_if(opt.faces == BuildOptions::FaceOption::none, "Grid::finalize: node_faces requires faces != none");
            const FacesCache& fc = *m_faces;
            const id_t nf = (id_t)fc.owner_cell.size();

            build_node_incidence(
                ic.node_faces, n_nodes,
                [&](auto push){
                    for (id_t fid = 0; fid < nf; ++fid) {
                        for (id_t nid : fc.face_nodes[fid]) push(nid, fid);
                    }
                }
            );
        }
    }

    // -----------------------------
    // 5) build EdgesCache (optional)
    // -----------------------------
    if (opt.build_edges) {
        m_edges.emplace();
        EdgesCache& ec = *m_edges;
        ec.clear();

        // key = (min<<32)|max (id_t expected <= 32-bit)
        auto edge_key = [&](id_t a, id_t b) -> std::uint64_t {
            id_t lo = std::min(a,b), hi = std::max(a,b);
            return (std::uint64_t(lo) << 32) | std::uint64_t(hi);
        };

        std::unordered_map<std::uint64_t, id_t> emap;
        emap.reserve((size_t)n_cells * 12);

        std::vector<std::vector<id_t>> edge_faces_rows;
        std::vector<std::vector<id_t>> edge_cells_rows;

        auto ensure_edge = [&](id_t a, id_t b) -> id_t {
            const std::uint64_t k = edge_key(a,b);
            auto it = emap.find(k);
            if (it != emap.end()) return it->second;

            const id_t eid = (id_t)ec.edges.size();
            emap.emplace(k, eid);
            ec.edges.push_back({ std::min(a,b), std::max(a,b) });
            edge_faces_rows.emplace_back();
            edge_cells_rows.emplace_back();
            return eid;
        };

        // Prefer edges from faces (if built), else derive from cells.
        if (opt.faces != BuildOptions::FaceOption::none) {
            const FacesCache& fc = *m_faces;
            const id_t nf = (id_t)fc.owner_cell.size();

            for (id_t fid = 0; fid < nf; ++fid) {
                std::vector<id_t> cyc = corners_from_face_nodes(std::vector<id_t>(fc.face_nodes[fid].begin(), fc.face_nodes[fid].end()));
                if (cyc.size() == 2) {
                    const id_t eid = ensure_edge(cyc[0], cyc[1]);
                    edge_faces_rows[eid].push_back(fid);
                    // add owner+neighbor (if any)
                    edge_cells_rows[eid].push_back(fc.owner_cell[fid]);
                    if (fc.neighbor_cell[fid] != invalid_id) edge_cells_rows[eid].push_back(fc.neighbor_cell[fid]);
                } else if (cyc.size() >= 3) {
                    for (size_t i = 0; i < cyc.size(); ++i) {
                        const id_t a = cyc[i];
                        const id_t b = cyc[(i+1) % cyc.size()];
                        const id_t eid = ensure_edge(a,b);
                        edge_faces_rows[eid].push_back(fid);
                        edge_cells_rows[eid].push_back(fc.owner_cell[fid]);
                        if (fc.neighbor_cell[fid] != invalid_id) edge_cells_rows[eid].push_back(fc.neighbor_cell[fid]);
                    }
                }
            }
        } else {
            // Faces not built: add edges from extracted tmp_cell_faces (still available because need_tmp_faces=true when build_edges)
            for (id_t cid = 0; cid < n_cells; ++cid) {
                for (const TempFace& tf : tmp_cell_faces[cid]) {
                    std::vector<id_t> cyc = corners_from_face_nodes(tf.nodes);
                    if (cyc.size() == 2) {
                        const id_t eid = ensure_edge(cyc[0], cyc[1]);
                        edge_cells_rows[eid].push_back(cid);
                    } else if (cyc.size() >= 3) {
                        for (size_t i = 0; i < cyc.size(); ++i) {
                            const id_t eid = ensure_edge(cyc[i], cyc[(i+1)%cyc.size()]);
                            edge_cells_rows[eid].push_back(cid);
                        }
                    }
                }
            }
        }

        // pack incidence
        build_csr(ec.edge_faces, edge_faces_rows);
        build_csr(ec.edge_cells, edge_cells_rows);

        ec.validate_or_throw();
    }

    // -----------------------------
    // 6) GeometryCache (optional)
    // -----------------------------
    if (opt.compute_face_geometry || opt.compute_cell_geometry) {
        m_geom.emplace();
        GeometryCache& gc = *m_geom;
        gc.clear();

        // face geometry
        if (opt.compute_face_geometry) {
            throw_if(opt.faces == BuildOptions::FaceOption::none, "Grid::finalize: face geometry requires faces != none");
            const FacesCache& fc = *m_faces;
            const id_t nf = (id_t)fc.owner_cell.size();

            gc.face_centroids.resize(nf);
            gc.face_normals.resize(nf);
            gc.face_areas.resize(nf);

            for (id_t fid = 0; fid < nf; ++fid) {
                std::vector<id_t> fnodes;
                fnodes.reserve(fc.face_nodes[fid].size());
                for (id_t nid : fc.face_nodes[fid]) fnodes.push_back(nid);

                const id_t owner = fc.owner_cell[fid];
                const Vector3d cc = cell_centroid_guess(owner);

                auto [n, a, fc0] = face_geom_from_nodes(fnodes, cc);
                gc.face_centroids[fid] = fc0;
                gc.face_normals[fid]   = n;
                gc.face_areas[fid]     = a;
            }
        }

        // cell geometry
        if (opt.compute_cell_geometry) {
            gc.cell_centroids.resize(n_cells, Vector3d(0,0,0));
            gc.cell_volumes.resize(n_cells, 0.0);

            // Compute via per-cell temp faces (robust for unique faces too).
            for (id_t cid = 0; cid < n_cells; ++cid) {
                const Vector3d cc_guess = cell_centroid_guess(cid);

                if (m_dim == 2) {
                    // Use boundary corner polygon depending on type.
                    const Cell& c = m_cells[cid];
                    std::vector<id_t> poly;

                    if (c.type == CellType::TRIANGLE && c.nodes.size() >= 3) {
                        poly = { c.nodes[0], c.nodes[1], c.nodes[2] };
                    } else if (c.type == CellType::QUAD && c.nodes.size() >= 4) {
                        poly = { c.nodes[0], c.nodes[1], c.nodes[2], c.nodes[3] };
                    } else if (c.type == CellType::POLYGON) {
                        poly.assign(c.nodes.begin(), c.nodes.end());
                    } else if (c.type == CellType::AMR2D && c.nodes.size() == 9) {
                        poly = { c.nodes[0], c.nodes[2], c.nodes[8], c.nodes[6] }; // corners
                    } else {
                        // fallback: use all nodes (may be wrong for strange ordering)
                        poly.assign(c.nodes.begin(), c.nodes.end());
                    }

                    // Triangulate fan around p0
                    const Vector3d p0 = node_pos(poly[0]);
                    double area = 0.0;
                    Vector3d cacc(0,0,0);

                    for (size_t i = 1; i + 1 < poly.size(); ++i) {
                        const Vector3d p1 = node_pos(poly[i]);
                        const Vector3d p2 = node_pos(poly[i+1]);
                        const Vector3d av = 0.5 * (p1 - p0).cross(p2 - p0);
                        const double a = av.norm(); // unsigned
                        area += a;
                        cacc += a * (p0 + p1 + p2) / 3.0;
                    }

                    if (area > 0.0) {
                        gc.cell_volumes[cid]   = area;          // "volume" = area in 2D
                        gc.cell_centroids[cid] = cacc / area;
                    } else {
                        gc.cell_volumes[cid]   = 0.0;
                        gc.cell_centroids[cid] = cc_guess;
                    }
                } else {
                    // 3D: signed volume/centroid from tetra fan to origin (0,0,0)
                    double V = 0.0;
                    Vector3d C(0,0,0);
                    const Vector3d O(0,0,0);

                    for (const TempFace& tf : tmp_cell_faces[cid]) {
                        std::vector<id_t> cyc = corners_from_face_nodes(tf.nodes);
                        if (cyc.size() < 3) continue;

                        const Vector3d a0 = node_pos(cyc[0]);
                        for (size_t i = 1; i + 1 < cyc.size(); ++i) {
                            const Vector3d a1 = node_pos(cyc[i]);
                            const Vector3d a2 = node_pos(cyc[i+1]);

                            // signed volume of tetra (O, a0, a1, a2)
                            const double v = a0.dot((a1).cross(a2)) / 6.0;
                            V += v;
                            const Vector3d tet_c = (O + a0 + a1 + a2) / 4.0;
                            C += v * tet_c;
                        }
                    }

                    if (std::abs(V) > 0.0) {
                        const double Vabs = std::abs(V);
                        gc.cell_volumes[cid]   = Vabs;
                        gc.cell_centroids[cid] = C / V; // keep sign consistent then ok
                    } else {
                        gc.cell_volumes[cid]   = 0.0;
                        gc.cell_centroids[cid] = cc_guess;
                    }
                }
            }
        }
    }

    // -----------------------------
    // 7) finalize state
    // -----------------------------
    d.clear();
    m_draft.reset();
    m_state = State::finalized;

    report.clear();
    if (bool valid = validate_finalized_basic(&report); !valid) {
        std::cerr << report << std::endl;
        throw std::runtime_error("Grid::finalize(): Bad finalized");
    }
}

} // namespace zephyr::geom