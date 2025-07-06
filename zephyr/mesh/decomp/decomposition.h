#pragma once

#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_prim.h>

namespace zephyr::utils { class Json; }

namespace zephyr::mesh::decomp {

using zephyr::geom::Box;

/// @brief Базовый класс декомпозиции
class Decomposition {
public:
    /// @brief Умный указатель на экземаляр класса
    using Ptr = std::shared_ptr<Decomposition>;
    using Ref = const std::shared_ptr<Decomposition>&;

    /// @brief Конструктор по умолчанию.
    /// Размер декомпозиции выбирается равным mpi::size()
    Decomposition();

    /// @brief Конструктор
    /// @param size Размер декомпозиции, не обязательно равен mpi::size(),
    /// для возможности тестирования декомпозиции на одном процессе.
    explicit Decomposition(int size);

    /// @brief Создать декомпозицию по конфигурации
    static Decomposition::Ptr create(const Box& domain, const utils::Json& config);

    /// @brief Деструктор по умолчанию
	virtual ~Decomposition() = default;

	/// @brief Размер декомпозиции
	int size() const { return m_size; }

    /// @brief Основная функция. Определение нового ранга ячейки.
    virtual int rank(EuCell& elem) const = 0;

	/// @brief Дисбаланс нагрузки
    static double imbalance(const std::vector<double>& ws);

	/// @brief Балансировка нагрузки
    virtual void balancing(const std::vector<double>& w) = 0;


protected:
    ///< Размер декомпозиции, не обязательно равен mpi::size(),
    /// для возможности тестирования декомпозиции на одном процессе.
    int m_size;
};

} // namespace zephyr::mesh::decomp