#pragma once

#include <zephyr/mesh/decomp/decomposition.h>
#include <zephyr/mesh/decomp/vdiagram.h>

namespace zephyr::mesh::decomp {

/// @brief Декомпозиция на основе диаграмм Вороного
class VD3 : public Decomposition {
public:
    /// @brief Умный указатель на экземаляр класса
    using Ptr = std::shared_ptr<VD3>;
    using Ref = const std::shared_ptr<VD3>&;

    /// @brief Конструктор со случайными генераторами
    /// @param size Размер декомпозиции, не обязательно равен mpi::size(),
    /// для возможности тестирования декомпозиции на одном процессе.
    VD3(const Box& domain, int size);

    /// @brief Создать указатель
    static VD3::Ptr create(const Box& domain, int size) {
        return std::make_shared<VD3>(domain, size);
    }
    
    /// @brief Основная функция. Определение нового ранга ячейки.
    int rank(const EuCell& elem) const final;

    /// @brief Балансировка нагрузки
    void balancing(const std::vector<double>& w) final;

protected:
    VDiagram m_diagram;
};

} // namespace zephyr::mesh::decomp