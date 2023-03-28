#include <zephyr/mesh/mesh.h>


namespace zephyr { namespace mesh {


void Mesh::initialize(Generator* gen) {
    gen->initialize(m_locals);
}



} // mesh
} // zephyr