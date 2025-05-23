#pragma once

#include <zephyr/mesh/mesh.h>

#include <zephyr/io/filter.h>
#include <zephyr/io/variables.h>


namespace zephyr::io {

using zephyr::mesh::AmrStorage;

/// @class VtuFile для записи хранилища в файл VTU для неструктурированных
/// сеток. Частицы также могут быть записаны.
class VtuFile {
public:

    /// @brief Все переменные класса имеют публичный доступ, запись файла
    /// можно осуществлять после настройки всех параметров

    std::string   filename;   ///< Полное имя файла
    Variables     variables;  ///< Список переменных на запись
    ComplexFilter filter;     ///< Функтор, возвращает true для нужных ячеек
    bool hex_only = true;     ///< Записывать как четырехугольники/шестигранники

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
    void save(mesh::SoaMesh &mesh);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::LaMesh &mesh);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(AmrStorage &cells);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::AmrCells &cells);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    /// @param nodes Массив уникальных вершин
    void save(AmrStorage &cells, const std::vector<geom::Vector3d>& nodes);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(CellStorage &cells, NodeStorage& nodes);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, AmrStorage &cells,
                     const Variables &variables, 
                     bool hex_only   = false,
                     bool polyhedral = false,
                     const Filter &filter = TrivialFilter());

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, mesh::AmrCells &locals,
                     const Variables &variables,
                     bool hex_only   = false,
                     bool polyhedral = false,
                     const Filter &filter = TrivialFilter());

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, AmrStorage &cells,
                     const std::vector<geom::Vector3d>& nodes,
                     const Variables &variables, 
                     bool hex_only   = false,
                     bool polyhedral = false);

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename, CellStorage &cells,
                     NodeStorage &nodes, const Variables &variables,
                     const Filter &filter = TrivialFilter());

    /// @brief Простое сохранение узлов
    static void save(const std::string &filename, NodeStorage &nodes);

};

} // namespace zephyr::io