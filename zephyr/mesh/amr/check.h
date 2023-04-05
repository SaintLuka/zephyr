/*

/// @file Файл содержит объявления функций, которые используются дебаге, а также
/// для проверки сетки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Выводит информацию о геометрии ячейки, списки вершин, граней,
/// соседей и прочее
void print_cell_info(Storage::Item cell);

/// @brief Выводит полную информацию о ячейке, а также о её окружении
void print_cell_info(Storage& locals, Storage& aliens, size_t ic);

/// @brief Выводит информацию о геометрии ячейки и её соседях в формате
/// скрипта для визуализации на python
void visualize_cell(Storage::Item cell);

/// @brief Проверяет наличие в Storage типов данных, необходимых для адаптации
/// @throw runtime_error, если один из необходимых типов отсутствует
void test_storage(Storage &cells);

/// @brief Проверяет базовую сетку на возможность адаптации, требуется наличие
/// в Storage всех необходимых типов, правильная геометрия ячеек (четырехугольники
/// в 2D или шестигранники типа куба в 3D), а также правильное расположение вершин
/// и граней в списках вершин и граней.
/// @return -1, если базовая сетка не подходящая
int check_base_mesh(Storage &locals, Storage &aliens, int rank = 0);

/// @brief Проверяет адаптированную сетку (связность ячеек, выполнение баланса и
/// так далее), функция используется при дебаге
/// @return -1, если сетка адаптирована неверно
int check_refined_mesh(Storage &locals, Storage &aliens, int rank = 0);

} // namespace amr
} // namespace mesh
} // namespace zephyr
*/