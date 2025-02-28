#include <iomanip>
#include <zephyr/utils/threads.h>

#include <zephyr/geom/sections.h>
#include <zephyr/geom/intersection.h>
#include <zephyr/geom/interface_recovery.h>

namespace zephyr::geom {

using namespace zephyr::utils;

inline double cross(const Vector3d& v1, const Vector3d& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

// Производная по теореме Гаусса--Остроградского
void InterfaceRecovery::compute_normal(EuCell &cell) const {
    double uc = cell(a);

    if (uc < 1.0e-12 || uc > 1.0 - 1.0e-12) {
        cell(n) = Vector3d::Zero();
        return;
    }

    cell(n) = Vector3d::Zero();
    for (auto &face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        double un = face.neib(a);

        Vector3d S = face.normal() * face.area();
        cell(n) -= face_fraction(uc, un) * S;
    }

    cell(n).normalize();
}

void InterfaceRecovery::compute_normals(EuMesh &mesh) const {
    for (auto cell: mesh) {
        compute_normal(cell);
    }
}

void InterfaceRecovery::find_section(EuCell &cell) const {
    if (cell(n).isZero()) {
        cell(p) = cell.center();
    } else {
        auto poly = cell.polygon();
        cell(p) = poly.find_section(cell(n), cell(a));
    }
}

void InterfaceRecovery::find_sections(EuMesh &mesh) const {
    for (auto cell: mesh) {
        find_section(cell);
    }
}

void InterfaceRecovery::adjust_normal(EuCell &cell) const {
    if (cell(n).isZero()) {
        return;
    }

    obj::plane plane{cell(p), cell(n)};

    std::vector <Vector3d> ints;
    for (auto face: cell.faces()) {
        obj::segment seg{face.vs(0), face.vs(1)};

        if (intersection2D::exist(plane, seg)) {
            if (face.neib(n).isZero()) {
                if (plane.n.dot(seg.tau()) * (face.neib(a) - 0.5) > 0.0) {
                    ints.emplace_back(seg.get(1.001));
                } else {
                    ints.emplace_back(seg.get(-0.001));
                }
            } else {
                obj::segment seg2{plane.p, face.neib(p)};
                Vector3d vi = intersection2D::find_fast(seg, seg2);
                ints.emplace_back(vi);
            }
        }
    }

    if (ints.size() > 1) {
        cell(n) = {ints[0].y() - ints[1].y(),
                       ints[1].x() - ints[0].x(), 0.0};

        if (cell(n).dot(plane.n) < 0.0) {
            cell(n) *= -1.0;
        }

        cell(n).normalize();
    }
}

void InterfaceRecovery::adjust_normals(EuMesh &mesh) const {
    for (auto cell: mesh) {
        adjust_normal(cell);
    }
}

void InterfaceRecovery::update(EuMesh &mesh, int smoothing) const {
    compute_normals(mesh);
    find_sections(mesh);
    for (int i = 0; i < smoothing; ++i) {
        adjust_normals(mesh);
        find_sections(mesh);
    }
}

AmrStorage InterfaceRecovery::body(EuMesh& mesh) const {
    int count = 0;
    for (auto cell: mesh) {
        if (cell(a) < 1.0e-12) {
            continue;
        }
        ++count;
    }

    AmrStorage cells(count);

    count = 0;
    for (auto cell: mesh) {
        if (cell(a) < 1.0e-12) {
            continue;
        }

        if (cell(a) < 1.0 - 1.0e-12) {
            if (cell(n).isZero()) {
                double d = 0.5 * std::sqrt(cell(a) * cell.volume());
                Quad quad = {
                        cell(p) + Vector3d{-d, -d, 0.0},
                        cell(p) + Vector3d{+d, -d, 0.0},
                        cell(p) + Vector3d{-d, +d, 0.0},
                        cell(p) + Vector3d{+d, +d, 0.0},
                };
                cells[count] = mesh::AmrCell(quad);
            }
            else {
                auto poly = cell.polygon();
                auto part = poly.clip(cell(p), cell(n));
                cells[count] = mesh::AmrCell(part);
            }
        }
        else {
            cells[count] = cell.geom();
        }
        ++count;
    }

    return cells;
}

} // namespace zephyr::geom