#pragma once

#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/io/filter.h>
#include <zephyr/io/variables.h>


namespace zephyr::io {

/// @class CsvFile для записи хранилища в файл .csv для неструктурированных
/// сеток.
class CsvFile {
    using AmrCells = zephyr::mesh::AmrCells;
public:

    /// @brief Все переменные класса имеют публичный доступ, запись файла
    /// можно осуществлять после настройки всех параметров

    std::string   filename;   ///< Полное имя файла
    int           precision;  ///< Точность записи
    Variables     variables;  ///< Список переменных на запись
    ComplexFilter filter;     ///< Функтор, возвращает true для нужных ячеек


    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    explicit CsvFile(const std::string &filename, int precision = 6,
                     const Variables &variables = {});

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(AmrCells &cells);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса write, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, AmrCells &cells,
                     int precision, const Variables &variables,
                     const Filter &filter = TrivialFilter());

};

} // namespace zephyr::io