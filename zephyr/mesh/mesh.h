#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/lagrange/la_mesh.h>

#include <zephyr/mesh/face.h>
#include <zephyr/mesh/cell.h>

namespace zephyr::mesh {

// В дальнейшем реализовать полиморфизм
using Mesh = EuMesh;

}