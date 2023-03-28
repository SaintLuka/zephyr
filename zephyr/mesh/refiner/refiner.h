#pragma once

#include <zephyr/data/storage.h>
#include <zephyr/mesh/refiner/distributor.h>

#ifdef ZEPHYR_ENABLE_MULTITHREADING
#include <zephyr/multithreading/thread-pool.h>
#endif

#ifdef ZEPHYR_ENABLE_MPI
#include <zephyr/network/mpi/network.h>
#include <zephyr/network/decomposition/base.h>
#endif

namespace zephyr { namespace mesh {

using zephyr::data::Storage;

#ifdef ZEPHYR_ENABLE_MULTITHREADING
using zephyr::multithreading::ThreadPool;
using zephyr::multithreading::dummy_pool;
#endif

#ifdef ZEPHYR_ENABLE_MPI
using ::zephyr::network::mpi::Network;
using Decomposition = zephyr::network::decomposition::Base;
#endif

/**
    \class
        \~russian Класс для использования адаптации, на данный момент содержит
        только статические функции, по сути можно заменить на namespace.
*/
class Refiner {
public:

    /**
	\brief
        \~russian Осуществляет инициализацию хранилища перед использованием
        функций адаптации, выполняется один раз после создания хранилища.
    */
    static void init(Storage &cells);

#ifdef ZEPHYR_ENABLE_MPI
    /**
	\brief
        \~russian Осуществляет инициализацию хранилища перед использованием
        функций адаптации, выполняется один раз после создания хранилища.
    */
    static void init(Network& network, Storage &cells);
#endif

    /**
	\brief Выводит полную информацию о ячейке
    */
    static void print_cell_info(Storage::iterator cell);

    /**
	\brief Выводит полную информаци о ячейке и её окружении
    */
    static void print_cell_info(Storage& locals, Storage& aliens, size_t ic);

    /**
	\brief Выводит полную информацию о ячейке в виде python
     скрипта для визуализации
    */
    static void visualize_cell(Storage::iterator cell);

    /**
	\brief
        \~russian Проверить базовую сетку после создания.
    \return
        \~russian -1, если сетка не подходит для адаптации.
    */
    static int check_base_mesh(Storage &cells);

#ifdef ZEPHYR_ENABLE_MPI
    /**
	\brief
        \~russian Проверить базовую сетку после создания.
    \return
        \~russian -1, если сетка не подходит для адаптации.
    */
    static int check_base_mesh(Decomposition& decomposition);
#endif

    /**
	\brief
        \~russian Проверить сетку после адаптации (для дебага).
    \return
        \~russian -1, если сетка имеет неверную структуру.
    */
    static int check_refined_mesh(Storage &cells);

#ifdef ZEPHYR_ENABLE_MPI
    /**
	\brief
        \~russian Проверить сетку после адаптации (для дебага).
    \return
        \~russian -1, если сетка имеет неверную структуру.
    */
    static int check_refined_mesh(Decomposition& decomposition);
#endif

    /**
	\brief
        \~russian Основная функция адаптации, меняет хранилище в соответствии
        с флагами amrData.flag ячеек.
    \param cells
        \~russian Хранилище ячеек.
    \param max_level
        \~russian Максимальный уровнь адаптации.
    \param op
        \~russian Оператор распределения физических данных.
    */
    static void refine(
            Storage &cells,
            unsigned int max_level,
            const DataDistributor& op
            if_multithreading(, ThreadPool& threads = dummy_pool)
    );

#ifdef ZEPHYR_ENABLE_MPI
    /**
    \brief
        \~russian Основная функция адаптации, меняет хранилище в соответствии
        с флагами amrData.flag ячеек.
    \param decomposition
        \~russian Декомпозиция хранилища.
    \param max_level
        \~russian Максимальный уровнь адаптации.
    \param op
        \~russian Оператор распределения физических данных.
    */
    static void refine(
            Decomposition& decomposition,
            unsigned int max_level,
            const DataDistributor& op
            if_multithreading(, ThreadPool& threads = dummy_pool)
    );

#endif
    /**
	\brief
        \~russian Основная функция адаптации, адаптация без распределения
        физических данных ячеек.
    \param cells
        \~russian Хранилище ячеек.
    \param max_level
        \~russian Максимальный уровнь адаптации.
    */
    static void refine(
            Storage &cells,
            unsigned int max_level
            if_multithreading(, ThreadPool& threads = dummy_pool)
    ) {
        refine(cells, max_level, { } if_multithreading(, threads));
    }

#ifdef ZEPHYR_ENABLE_MPI
    /**
    \brief
        \~russian Основная функция адаптации, адаптация без распределения
        физических данных ячеек
    \param decomposition
        \~russian Декомпозиция хранилища.
    \param max_level
        \~russian Максимальный уровнь адаптации.
    */
    static void refine(
            Decomposition& decomposition,
            unsigned int max_level
            if_multithreading(, ThreadPool& threads = dummy_pool)
    ) {
        refine(decomposition, max_level, { } if_multithreading(, threads));
    }
#endif
};

} // namespace mesh
} // namespace zephyr