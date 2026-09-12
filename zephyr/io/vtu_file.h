#pragma once

#include <zephyr/io/variables.h>

// Forward declaration
namespace zephyr::mesh {
class RawCells;
class RawNodes;
class Mesh;
}

namespace zephyr::io {

/// @brief Опции записи
struct VtuOptions {
    /// @brief Для адаптивных сеток: записывать ячейки как полигоны? В обратном
    /// случае ячейки пишутся как простые четырехугольники/шестигранники,
    /// то есть с "висящими" узлами.
    ///
    /// Для трёхмерных сеток: интерпретировать ячейки как многогранники общего вида.
    /// Сохраняет отдельно грани. Работает дольше, файлы весят больше.
    /// Необходимо использовать при записи сетки с многогранниками.
    bool polyhedral = false;

    /// @brief Сохранять уникальные вершины? Если сетка уже содержит уникальные
    /// вершины, то использует их. Если сетка эйлерова, то уникальные вершины
    /// генерируются, в этом случае запись займет больше времени.
    bool unique_nodes = false;
};

/// @brief Запись неструктурированной сетки в VTU-файл.
class VtuFile {
public:
    // Переменные класса имеют публичный доступ, сохранение можно выполнить
    // после настройки всех параметров

    std::string filename;   ///< Полное имя файла
    Variables   variables;  ///< Список переменных на запись
    VtuOptions  options;    ///< Опции записи

    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    explicit VtuFile(std::string_view filename,
                     const Variables& variables = {},
                     const VtuOptions& options = {});

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::Mesh& mesh) const;

    /// @brief Базовая функция записи в файл. До вызова функции должен быть
    /// создан экземпляр класса и настроены опции записи.
    void save(mesh::RawCells& cells) const;

    /// @brief Статическая функция записи в файл. Полный аналог функции-члена
    /// класса save, но вызывается без экземпляра класса, все параметры записи
    /// передаются непосредственно как аргументы функции.
    static void save(std::string_view filename,
                     mesh::Mesh& mesh,
                     const Variables& variables = {},
                     const VtuOptions& options = {});

    static void save(std::string_view filename,
                     mesh::RawCells& locals,
                     const Variables& variables = {},
                     const VtuOptions& options = {});

};

} // namespace zephyr::io