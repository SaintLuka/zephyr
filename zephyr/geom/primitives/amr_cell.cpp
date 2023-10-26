#include <zephyr/geom/primitives/amr_cell.h>

namespace zephyr::geom {

AmrCell::AmrCell(const Quad& quad)
    : Element(0, 0), dim(2), linear(true),
    vertices(quad), faces(2),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double area = quad.area();

    size = std::sqrt(area);
    center = quad.centroid(area);

    for (int i = 0; i < 4; ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqQuad& quad)
    : AmrCell(quad.reduce()) { }

AmrCell::AmrCell(const Cube& cube)
    : Element(0, 0), dim(3), linear(true),
    vertices(cube), faces(3),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double volume = cube.volume();
    size = std::cbrt(volume);
    center = cube.centroid(volume);

    for (int i = 0; i < 6; ++i) {
        Quad vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]],
                vertices[faces[i].vertices[2]],
                vertices[faces[i].vertices[3]]
        };

        faces[i].area     = vs.area();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqCube& cube)
    : AmrCell(cube.reduce()) {

}

double AmrCell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

} // namespace zephyr::geom