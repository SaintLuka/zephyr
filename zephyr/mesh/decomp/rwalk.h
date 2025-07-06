#pragma once

#include <zephyr/mesh/decomp/decomposition.h>
#include <zephyr/mesh/decomp/vdiagram.h>

namespace zephyr::mesh::decomp {

/// @brief RWalk от Random Walk, блуждающие ячейки Вороного.
class RWalk : public Decomposition {
public:
    /// @brief Умный указатель на экземаляр класса
    using Ptr = std::shared_ptr<RWalk>;
    using Ref = const std::shared_ptr<RWalk>&;

    /// @brief Конструктор со случайными генераторами
    /// @param size Размер декомпозиции, не обязательно равен mpi::size(),
    /// для возможности тестирования декомпозиции на одном процессе.
    RWalk(const Box& domain, int size);

    /// @brief Создать указатель
    static RWalk::Ptr create(const Box& domain, int size) {
        return std::make_shared<RWalk>(domain, size);
    }

    /// @brief Основная функция. Определение нового ранга ячейки.
    int rank(EuCell& elem) const final;

    /// @brief Балансировка нагрузки (ничего не балансирует,
    /// просто смещает генераторы ячеек)
    void balancing(const std::vector<double>& w) final;

protected:
    double   m_step;
    Box      m_domain;
    VDiagram m_diagram;
};

} // namespace zephyr::mesh::decomp