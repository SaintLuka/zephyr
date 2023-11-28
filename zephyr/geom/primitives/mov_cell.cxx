#include <zephyr/geom/primitives/mov_cell.h>

namespace zephyr::geom {

MovCell::MovCell(const Polygon& poly)
        : Element(0, 0), dim(2),
          faces(CellType::POLYGON, poly.size()) {

    Vector3d C = poly.center();
    double area = poly.area(C);

    size = std::sqrt(area);
    center = poly.centroid(area);

    for (int i = 0; i < poly.size(); ++i) {
        Line vs = {
                poly[faces[i].vertices[0]],
                poly[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

double MovCell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

} // namespace zephyr::geom