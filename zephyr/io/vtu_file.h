#pragma once

#include <zephyr/io/variables.h>

// Forward declaration
namespace zephyr::mesh {
class AmrCells;
class EuMesh;
}

namespace zephyr::io {

/// @class VtuFile для записи хранилища в файл VTU для неструктурированных
/// сеток. Частицы также могут быть записаны.
class VtuFile {
public:

    /// @brief Все переменные класса имеют публичный доступ, запись файла
    /// можно осуществлять после настройки всех параметров

    std::string filename;   ///< Полное имя файла
    Variables   variables;  ///< Список переменных на запись
    bool hex_only = true;   ///< Записывать как четырехугольники/шестигранники

    /// @brief Интерпретировать трёхмерные ячейки как многогранники общего вида.
    /// Сохраняет отдельно грани. Работет долго, включать при необходимости.
    bool polyhedral = false;
    
    /// @brief Сохранять уникальные вершины, актуально для EuMesh, если сетка
    /// часто перестраивается, то функция выполняется достаточно долго.
    bool unique_nodes = false;


    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    explicit VtuFile(const std::string &filename,
                     const Variables &variables = {},
                     bool hex_only   = true,
                     bool polyhedral = false);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::EuMesh &mesh);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::AmrCells &cells);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename,
                     mesh::AmrCells &locals,
                     const Variables &variables,
                     bool hex_only   = false,
                     bool polyhedral = false);

};

} // namespace zephyr::io