#pragma once

#include <zephyr/io/variables.h>

// Forward declaration
namespace zephyr::mesh {
class AmrCells;
class AmrNodes;
class EuMesh;
}

namespace zephyr::io {

/// @brief Запись неструктурированной сетки в VTU-файл.
class VtuFile {
public:
    // Переменные класса имеют публичный доступ, сохранение можно выполнить
    // после настройки всех параметров

    std::string filename;   ///< Полное имя файла
    Variables   variables;  ///< Список переменных на запись

    /// @brief Для адаптивных сеток: записывать ячейки как полигоны? В обратном
    /// случае ячейки пишутся как простые четырехугольники/шестигранники,
    /// то есть с "висящими" узлами.
    ///
    /// Для трёхмерных сеток: интерпретировать ячейки как многогранники общего вида.
    /// Сохраняет отдельно грани. Работает долго, включать при необходимости.
    /// Необходимо использовать при записи сетки с многогранниками.
    bool polyhedral = false;
    
    /// @brief Сохранять уникальные вершины, актуально для EuMesh. Если сетка
    /// часто перестраивается, то функция выполняется достаточно долго.
    bool unique_nodes = false;


    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    explicit VtuFile(const std::string &filename,
                     const Variables &variables = {},
                     bool polyhedral   = false,
                     bool unique_nodes = false);

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::EuMesh &mesh) const;

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::AmrCells &cells) const;

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::AmrCells &cells, const mesh::AmrNodes& nodes) const;

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(const std::string &filename,
                     mesh::EuMesh &mesh,
                     const Variables &variables = {},
                     bool polyhedral = false,
                     bool unique_nodes = false);

    static void save(const std::string &filename,
                     mesh::AmrCells &locals,
                     const Variables &variables = {},
                     bool polyhedral = false);

    static void save(const std::string &filename,
                     mesh::AmrCells &locals,
                     const mesh::AmrNodes &nodes,
                     const Variables &variables = {},
                     bool polyhedral = false);

};

} // namespace zephyr::io