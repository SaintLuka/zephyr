#include <zephyr/mesh/cell.h>
#include <zephyr/mesh/range.h>

namespace zephyr { namespace mesh {

Range::Range(Storage &locals, Storage &aliens, int from, int to, int id)
    : m_id(id), m_locals(locals), m_aliens(aliens) {

    m_begin = std::min(from, locals.size());
    m_end = std::min(to, m_locals.size());

    // Пустой диапозон
    if (m_begin >= m_end) {
        m_begin = 0;
        m_end = 0;
    }
}

int Range::size() const {
    return m_end - m_begin;
}

int Range::id() const {
    return m_id;
}

ICell Range::begin() const {
    return { m_locals, m_aliens, m_begin };
}

ICell Range::end() const {
    return { m_locals, m_aliens, m_end };
}

} // namespace mesh
} // namespace zephyr