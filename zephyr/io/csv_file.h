#pragma once

#include <zephyr/mesh/storage.h>
#include <zephyr/io/filter.h>
#include <zephyr/io/variables.h>


namespace zephyr { namespace io {

using zephyr::mesh::Storage;

/// @class CsvFile для записи хранилища в файл .csv для неструктурированных
/// сеток. Частицы также могут быть записаны.
class CsvFile {
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
    void save(Storage &cells);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса write, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, Storage &cells,
                     int precision, const Variables &variables,
                     const Filter &filter = TrivialFilter());

};

} // io
} // zephyr