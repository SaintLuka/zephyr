#pragma once

#include <zephyr/io/variable.h>

namespace zephyr::io {

/// @brief Список переменных для записи в VTU файл
/// Данный класс содержит массив дескрипторов полей для записи, подробнее
/// о дескрипторах см. в классе FieldDescriptor.
class Variables {
public:

    /// @brief Пустой список переменных
    Variables() = default;

    /// @brief Создать список по одному имени VariablesList = "uid"
    explicit Variables(const char *name);

    /// @brief Создать список по одному имени VariablesList = "uid"
    explicit Variables(const std::string &name);

    /// @brief Создать список по набору имен переменных
    /// Пример. VariablesList list = {"base_id", "coords", "level"}
    Variables(std::initializer_list<const char *> names);

    /// @brief Создать список по набору имен переменных
    /// Пример. VariablesList list = {"base_id", "coords", "level"}
    Variables(std::initializer_list<std::string> names);

    /// @brief Создать список по набору имен переменных
    /// Пример. VariablesList list = {"base_id", "coords", "level"}
    Variables(const std::vector<const char *> &names);

    /// @brief Создать список по набору имен переменных
    /// Пример. VariablesList list = {"base_id", "coords", "level"}
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
    /// переменные. К примеру, list.append("base_id"); list.append("level").
    void append(const char *name);

    /// @brief Аналогично функции append(const char* )
    void append(const std::string &name);

    /// @brief Добавть в список переменную с полным описанием
    /// @param name Имя переменной
    /// @param n_components Размер вектора для хранения переменной
    /// @param func Функция записи переменной по заданному указателю
    /// @tparam T Тип шаблона важен для правильной дедукции VtkType
    /// Следующий код позволяет добавить векторную переменную "momentum",
    /// которая позволяет записывать две компоненты импульса в формате
    /// Float32. В целом, использование чисел с одинарной точностью
    /// позволяет сократить размер выходного файла.
    /// @code
    /// VariableList list;
    /// list.append("momentum", 2,
    ///     WriteFunction<float>([](AmrStorage::Item& cell, float* out) {
    ///         out[0] = static_cast<float>(cell.mass * cell.velocity.x);
    ///         out[1] = static_cast<float>(cell.mass * cell.velocity.y);
    ///     }));
    /// @endcode

    template<class T>
    void append(const char *name, int n_components, const WriteSoaCell<T> &func) {
        m_list.emplace_back(name, n_components, func);
    }

    /// @brief Аналогично функции append(const char*, ...)
    template<class T>
    void append(const std::string &name, int n_components, const WriteSoaCell<T> &func) {
        m_list.emplace_back(name, n_components, func);
    }

    /// @brief Аналогично функции append(const char*, ... )
    template<class T>
    void append(const char *name, const WriteSoaCell<T> &func) {
        m_list.emplace_back(name, 1, func);
    }

    /// @brief Аналогично функции append(const char*, ...)
    template<class T>
    void append(const std::string &name, const WriteSoaCell<T> &func) {
        m_list.emplace_back(name, 1, func);
    }

    /// @brief Упрощенный вариант для добавления double полей
    void append(const std::string& name, std::function<double(mesh::EuCell&)> f) {
        m_list.emplace_back(name, 1, WriteSoaCell<double>(
                [f](mesh::EuCell& cell, double *out) {
                    out[0] = f(cell);
                }));
    }

    /// @brief Упрощенный синтаксис для добавления полей типа double
    /// @param name Название переменной
    /// @param func Некоторая функция, которая для произвольной ячейки выдает
    /// значение переменной
    void operator+=(std::pair<std::string, std::function<double(mesh::EuCell&)>> p) {
        append(p.first, p.second);
    }

    /// @brief Упрощенный синтаксис для добавления Storable полей
    /// @param name Название переменной
    /// @param func Некоторая функция, которая для произвольной ячейки выдает
    /// значение переменной
    template <typename T>
    void append(std::string name, const mesh::Storable<T>& p) {
        if (VtkType::get<T>().is_undefined()) {
            if constexpr (std::is_same_v<T, geom::Vector3d>) {
                append(name, 3, WriteSoaCell<double>(
                        [p](mesh::EuCell &cell, double *out) {
                            out[0] = cell(p).x();
                            out[1] = cell(p).y();
                            out[2] = cell(p).z();
                        }));
            } else {
                throw std::runtime_error("Can't add as VTK type");
            }
        }
        append(name, 1, WriteSoaCell<T>(
                [p](mesh::EuCell &cell, T *out) {
                    out[0] = cell(p);
                }));
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
    std::vector<Variable> m_list;
};

} // namespace zephyr::io

