#pragma once

#include <zephyr/mesh/generator/base.h>
#include <zephyr/mesh/storage.h>

#include <zephyr/mesh/range.h>


namespace zephyr { namespace mesh {


class Mesh {
public:
    using Generator = generator::Base;


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

    void initialize(Generator* gen);

private:
    Storage m_locals;
    Storage m_aliens;
};


} // mesh
} // zephyr