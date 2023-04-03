#pragma once

#include <zephyr/mesh/generator/generator.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/distributor.h>

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


    /// @brief Проверить базовую сетку после создания.
    /// @return -1, если сетка не подходит для адаптации.
    int check_base();

    /// @brief Проверить сетку после адаптации (для дебага).
    /// @return -1, если сетка имеет неверную структуру.
    int check_refined();

    /// @brief Основная функция адаптации, меняет хранилище в соответствии
    /// с флагами amr.flag ячеек.
    void refine();

private:
    void initialize(Generator* gen);

    /// @brief Осуществляет инициализацию хранилища перед использованием
    /// функций адаптации, выполняется один раз после создания хранилища.
    void init_amr();

    int max_level;
    Distributor distributor;

    Storage m_locals;
    Storage m_aliens;
};


} // mesh
} // zephyr