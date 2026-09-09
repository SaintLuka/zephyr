#pragma once

#include <zephyr/io/variable.h>
#include <zephyr/geom/vector.h>

namespace zephyr::io {

/// @brief Список переменных для записи в VTU файл
/// Данный класс содержит массив дескрипторов полей для записи, подробнее
/// о дескрипторах см. в классе FieldDescriptor.
class Variables {
public:

    /// @brief Пустой список переменных
    Variables() = default;

    /// @brief Создать список по одному имени Variables = "level"
    explicit Variables(std::string_view name);

    /// @brief Создать список по набору имен переменных
    /// Пример. Variables list = {"center", "level"}
    Variables(std::initializer_list<const char *> names);

    /// @brief Создать список по набору имен переменных
    /// Пример. Variables list = {"center", "level"}
    Variables(std::initializer_list<std::string> names);

    /// @brief Создать список по набору имен переменных
    /// Пример. Variables list = {"center", "level"}
    Variables(const std::vector<const char *> &names);

    /// @brief Создать список по набору имен переменных
    /// Пример. Variables list = {"center", "level"}
    Variables(const std::vector<std::string> &names);

    /// @brief Добавляет в список набор переменных по именам
    void append(std::initializer_list<const char *> names);

    /// @brief Добавляет в список набор переменных по именам
    void append(std::initializer_list<std::string> names);

    /// @brief Добавляет в список набор переменных по именам
    void append(const std::vector<const char *> &names);

    /// @brief Добавляет в список набор переменных по именам
    void append(const std::vector<std::string> &names);

    /// @brief Добавляет в список набор переменных, скопированных из другого
    /// списка переменных
    void append(const Variables &names);

    /// @details Добавить в список переменную по имени.
    /// Данным образом в список можно добавить только некоторые предопределенные
    /// переменные. К примеру, list.append("center"); list.append("level").
    void append(std::string_view name);

    /// @brief Добавить в список переменную с полным описанием
    /// @param name Имя переменной
    /// @param n_components Размер вектора для хранения переменной
    /// @param func Функция записи переменной по заданному указателю
    /// @tparam T Тип шаблона важен для правильной дедукции VtkType
    ///
    /// Следующий код добавляет векторную переменную "momentum", которая
    /// позволяет записывать две компоненты импульса в формате Float32.
    /// @code
    /// Variables vars;
    /// vars.append("momentum", 2,
    ///     WriteFunction<float>([](EuCell& cell, float* out) {
    ///         out[0] = static_cast<float>(cell.mass * cell.velocity.x);
    ///         out[1] = static_cast<float>(cell.mass * cell.velocity.y);
    ///     }));
    /// @endcode
    template<class T>
    void append(std::string_view name, int n_components, const WriteCell<T> &func) {
        list_.emplace_back(name, n_components, func);
    }

    template<class T>
    void append(std::string_view name, int n_components, const WriteNode<T> &func) {
        list_.emplace_back(name, n_components, func);
    }

    /// @brief Аналогично функции append(const char*, ... )
    template<class T>
    void append(std::string_view name, const WriteCell<T> &func) {
        list_.emplace_back(name, 1, func);
    }

    template<class T>
    void append(std::string_view name, const WriteNode<T> &func) {
        list_.emplace_back(name, 1, func);
    }

    /// @brief Упрощенный вариант для добавления скалярных полей ячеек
    template <typename T = double>
    std::enable_if_t<std::is_arithmetic_v<T>, void>
    append(std::string_view name, std::function<T(mesh::EuCell&)> func) {
        list_.emplace_back(name, 1, WriteCell<T>(
                [func](mesh::EuCell& cell, T *out) {
                    out[0] = func(cell);
                }));
    }

    /// @brief Упрощенный вариант для добавления скалярных полей узлов
    template <typename T = double>
    std::enable_if_t<std::is_arithmetic_v<T>, void>
    append(std::string_view name, std::function<T(mesh::EuNode&)> func) {
        list_.emplace_back(name, 1, WriteNode<T>(
                [func](mesh::EuNode& node, T *out) {
                    out[0] = func(node);
                }));
    }

    /// @brief Упрощенный вариант для добавления векторных полей ячеек
    template <typename T>
    std::enable_if_t<std::is_same_v<T, geom::Vector3d>, void>
    append(std::string_view name, std::function<geom::Vector3d(mesh::EuCell&)> func) {
        list_.emplace_back(name, 3, WriteCell<double>(
                [func](mesh::EuCell& cell, double *out) {
                    *reinterpret_cast<geom::Vector3d*>(out) = func(cell);
                }));
    }

    /// @brief Упрощенный вариант для добавления векторных полей узлов
    template <typename T>
    std::enable_if_t<std::is_same_v<T, geom::Vector3d>, void>
    append(std::string_view name, std::function<geom::Vector3d(mesh::EuNode&)> func) {
        list_.emplace_back(name, 3, WriteNode<double>(
                [func](mesh::EuNode& node, double *out) {
                    *reinterpret_cast<geom::Vector3d*>(out) = func(node);
                }));
    }

    /// @brief Упрощенный синтаксис для добавления полей типа double
    /// Variables vars;
    /// vars += {"rho", [](EuCell& cell) -> double { ... } };
    void operator+=(std::pair<std::string_view, std::function<double(mesh::EuCell&)>> p) {
        append(p.first, p.second);
    }

    /// @brief Упрощенный синтаксис для добавления полей типа double
    /// Variables vars;
    /// vars += {"rho", [](EuNode& node) -> double { ... } };
    void operator+=(std::pair<std::string_view, std::function<double(mesh::EuNode&)>> p) {
        append(p.first, p.second);
    }

    /// @brief Добавить существующий список к текущему
    void operator+=(const Variables& other) {
        for (const auto& var: other.list_) { list_.emplace_back(var); }
    }

    /// @brief Упрощенный синтаксис для добавления Storable полей типа int/double/Vector3d
    /// Storable<T> rho;
    /// Variables vars;
    /// vars += {"rho", rho};
    template <typename T>
    void add_cell_data(std::string_view name, mesh::Storable<T> p) {
        if (VtkType::get<T>().is_undefined()) {
            if constexpr (std::is_same_v<T, geom::Vector3d>) {
                append(name, 3, WriteCell<double>(
                        [p](mesh::EuCell &cell, double *out) {
                            out[0] = cell[p].x();
                            out[1] = cell[p].y();
                            out[2] = cell[p].z();
                        }));
            } else {
                throw std::runtime_error("Can't add as VTK type");
            }
        }
        else {
            append(name, 1, WriteCell<T>(
                [p](mesh::EuCell &cell, T *out) {
                    out[0] = cell[p];
                }));
        }
    }

    /// @brief Упрощенный синтаксис для добавления Storable полей типа int/double/Vector3d
    /// Storable<T> rho;
    /// Variables vars;
    /// vars += {"rho", rho};
    template <typename T>
    void add_node_data(std::string_view name, mesh::Storable<T> p) {
        if (VtkType::get<T>().is_undefined()) {
            if constexpr (std::is_same_v<T, geom::Vector3d>) {
                append(name, 3, WriteNode<double>(
                        [p](mesh::EuNode &node, double *out) {
                            out[0] = node[p].x();
                            out[1] = node[p].y();
                            out[2] = node[p].z();
                        }));
            } else {
                throw std::runtime_error("Can't add as VTK type");
            }
        }
        else {
            append(name, 1, WriteNode<T>(
                [p](mesh::EuNode &node, T *out) {
                    out[0] = node[p];
                }));
        }
    }

    /// @brief Очистить список переменных
    void reset();

    /// @brief Доступ к переменной в списке по индексу
    const Variable &operator[](int i) const;
    
    /// @brief Количество переменных
    size_t size() const;

    /// @brief Доступ ко всему списку переменных
    const std::vector<Variable> &list() const;

private:
    std::vector<Variable> list_;  ///< Список переменных
};

} // namespace zephyr::io