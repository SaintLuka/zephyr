#pragma once

#include <zephyr/mesh/generator/generator.h>
#include <zephyr/mesh/storage.h>

#include <zephyr/mesh/range.h>


namespace zephyr { namespace mesh {

using namespace generator;

class Mesh {
public:

    template <class T>
    Mesh(const T& val, Generator* gen)
            : m_locals(val), m_aliens(val) {
        initialize(gen);
    }


    Range cells() {
        return { m_locals, m_aliens, 0, m_locals.size() };
    }

    operator Storage&() { return m_locals; }

    Storage& locals() { return m_locals; }

    Storage& aliens() { return m_aliens; }


private:
    void initialize(Generator* gen);

    Storage m_locals;
    Storage m_aliens;
};


} // mesh
} // zephyr