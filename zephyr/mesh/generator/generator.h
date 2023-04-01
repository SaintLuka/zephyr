#pragma once

#include <memory>
#include <vector>

#include <zephyr/mesh/storage.h>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/face.h>

#ifdef ZEPHYR_ENABLE_YAML
#include <yaml-cpp/yaml.h>
#endif

namespace zephyr { namespace mesh { namespace generator {

using zephyr::geom::Vector3d;
using zephyr::geom::Box;
using zephyr::geom::FaceFlag;


/// @class Базовый класс для сеточных генераторов. Содержит набор
/// виртуальных функций для создания Storage с сеткой.
class Generator {
public:
    using Ptr = std::unique_ptr<Generator>;

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Создает один из генераторов
    static std::unique_ptr<Generator> create(YAML::Node config);
#endif

    /// @brief Базовый конструктор только с именем типа сетки
    explicit Generator(const std::string &type);

    /// @brief Тип сеточного генератора
    const std::string &type() const;

    /// @brief Количество элементов сетки
    virtual int size() const;

    /// @brief Ограничивающий объем
    virtual Box bbox() const;

    /// @brief Инициализация переданной сетки
    virtual void initialize(Storage& cells) = 0;

protected:

    /// @struct Прострая структура для описания части сетки
    struct Part {
    public:
        int rank;
        int n_cells;
        int n_parts;
        int from, to;

        /// @param n_cells Суммарное число ячеек
        /// @param n_parts Количество частей сетки
        /// @param rank Номер части (обычно rank)
        explicit Part(int n_cells, int n_parts = 1, int rank = 0);

        /// @brief Размер части сетки
        inline int size() const;
    };

    /// @brief Проверить размеры сетки перед созданием
    virtual void check_size() const;

    /// @brief Проверить параметры сетки перед созданием
    virtual void check_params() const;

    /// @brief Инициализация части переданной сетки
    virtual void initialize(Storage& storage, Part part);

    std::string m_type; ///< Имя сеточного генератора
    int m_size;         ///< Количество элементов сетки
};

} // generator
} // mesh
} // zephyr
