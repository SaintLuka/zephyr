#pragma once

#include <zephyr/data/storage.h>
#ifdef ZEPHYR_ENABLE_MULTITHREADING
#include <zephyr/multithreading/thread-pool.h>
#endif
#include <zephyr/network/decomposition/base.h>
#include <zephyr/network/decomposition/ORB/blocks.h>

namespace zephyr { namespace network { namespace decomposition {

#ifdef ZEPHYR_ENABLE_MULTITHREADING
using ::zephyr::multithreading::ThreadPool;
using ::zephyr::multithreading::dummy_pool;
#endif


class ORB : public Base {
public:

    /// @brief Значения по умолчанию для параметров балансировки
    struct defaults {
        /// @brief Использовать метод Ньютона
        constexpr static const bool newton = true;

        /// @brief Скрость смещения генераторов (0, 0.5)
        constexpr static const double mobility = 0.1;
    };

    ORB(Network&  network,
        Storage&  elements,
        Domain&   domain,
        Measurer& measurer,
        if_multithreading(ThreadPool& pool = dummy_pool,)
        bool   newton   = defaults::newton,
        double mobility = defaults::mobility
    );

#ifdef ZEPHYR_ENABLE_YAML
    ORB(Network&  network,
        Storage&  elements,
        Domain&   domain,
        Measurer& measurer,
        const YAML::Node& config
        if_multithreading(, ThreadPool& pool = dummy_pool)
    );
#endif

    // set функции

    void use_newton(bool val = true);

    void set_mobility(double val);

    // get функции

    bool newton() const;

    double mobility() const;

    Vector3d center() override;

    void balancing(const std::vector<double>& w) override;

protected:

    void collect_local_info() override;

    void collect_aliens_info() override;

    /// @brief Важная функция. Определить ранг процесса, которому принадлежит
    /// произвольная точка v. На практике точка v обычно является положением
    /// частицы или центром расчетной ячейки.
    /// Функция используется для перераспределения элементов между процессами.
    /// @param v Произвольная точка всей области
    int rank(const Vector3d& v) const final;

    /// @brief Важная функция. Определить, принадлежит ли некоторая точка v
    /// данного (!) процесса в некоторой окрестности процесса с рангом
    /// neib_rank. Функция используется для формирования списков aliens.
    /// В качестве радиуса поиска обычно выбирается максимальный радиус
    /// поиска для элементов с данного процесса и с соседнего (neib_rank)
    /// @param v Точка из подобласти данного процесса
    /// @param neib_rank Ранг процесса, в окрестности которого ищется точка
    /// @return true если точка принадлежит окрестности
    bool is_near(const Vector3d& v, int neib_rank) const final;

private:

    /// @brief Использовать метод Ньютона
    bool m_newton;

    /// @brief Скорость сдивга (0, 0.5)
    double m_mobility;

    /// @brief Блочная структура
    std::unique_ptr<Blocks> blocks;

    /// @brief Радиус поиска соседей
    std::vector<double> r_search;
};

} // decomposition
} // network
} // zephyr
