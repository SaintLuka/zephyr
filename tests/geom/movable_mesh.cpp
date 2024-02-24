#include <iostream>
#include <vector>
#include <complex>

#include <zephyr/mesh/mesh.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/io/pvd_file.h>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/cubic.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/block_structured.h>


using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::mesh::generator;


class _U_ {
public:
    double x;
    double y;
    double z;

    _U_() = default;
};

class _V_ {
public:
    Vector3d v;

    _V_() = default;
};


static const _U_ U{};
static const _V_ V{};

int main() {
    Rectangle rect(0.0, 1.0, 0.0, 1.0);
    rect.set_nx(200);

    LaMesh mesh(U, V, &rect);

    for (auto& cell: mesh.locals()) {
        cell(U).x = std::cos(3.0 / (cell.center.norm() + 0.1));
    }
    for (auto& node: mesh.nodes()) {
        node(U).x = std::cos(3.0 / (node.coords.norm() + 0.1));

        std::complex<double> z = {node.coords.x(), node.coords.y()};
        z *= z;

        Vector3d v = {z.real(), z.imag(), 0.0};

        node.shift.z() = v.norm(); //v - node.coords;
    }

    for (auto& node: mesh.nodes()) {
        node.move();
    }


    PvdFile pvd("rect", "out");

    pvd.variables += {"u", [](CellStorage::Item& cell) -> double { return cell(U).x; }};
    pvd.variables += {"v", [](NodeStorage::Item& node) -> double { return node(V).v.x(); }};

    pvd.save(mesh, 0.0);

    return 0;
}