#pragma once

#include <zephyr/mesh/storage.h>
#include <zephyr/io/filter.h>
#include <zephyr/io/variables.h>


namespace zephyr { namespace io {

using zephyr::mesh::Storage;

/// @class VtuFile для записи хранилища в файл VTU для неструктурированных
/// сеток. Частицы также могут быть записаны.
class VtuFile {
public:

    /// @brief Все переменные класса имеют публичный доступ, запись файла
    /// можно осуществлять после настройки всех параметров

    std::string   filename;  ///< Полное имя файла
    Variables variables; ///< Список переменных на запись
    ComplexFilter filter;    ///< Функтор, возвращает true для нужных ячеек
    bool          hex_only;  ///< Записывать как четырехугольники/шестигранники


    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    explicit VtuFile(const std::string &filename,
                     const Variables &variables = {},
                     bool hex_only = false);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(Storage &cells);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса write, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, Storage &cells,
                     const Variables &variables, bool hex_only = false,
                     const Filter &filter = TrivialFilter());

};

} // io
} // zephyr