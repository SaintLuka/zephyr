#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include <zephyr/io/vtk_type.h>

namespace zephyr::io {

/// @brief Запись структурированной сетки в VTR-файл.
class VtrFile {
public:
    /// @brief Конструктор класса, получает полный набор параметров.
    /// Также параметры могут быть заданы/изменены напрямую после создания
    /// экземпляра класса.
    VtrFile(const std::string &filename,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z = {});

    /// @brief Добавить поле данных для двумерной сетки
    template <typename T>
    std::enable_if_t<std::is_arithmetic_v<T>, void>
    point_data(const std::string& name, std::vector<std::vector<T>>& data);

    /// @brief Добавить поле данных для трёхмерной сетки
    template <typename T>
    std::enable_if_t<std::is_arithmetic_v<T>, void>
    point_data(const std::string& name, std::vector<std::vector<std::vector<T>>>& data);

    /// @brief Функция записи в файл
    void save() const;

protected:

    /// @brief Поле данных
    struct Field {
        VtkType type;
        std::vector<std::byte> data;
    };

    /// @brief Полное имя файла
    std::string m_filename;

    /// @brief Размерность сетки
    int m_dimension;

    /// @brief Массивы координат узлов
    std::vector<double> m_x, m_y, m_z;

    /// @brief Поля данных для сохранения
    std::unordered_map<std::string, Field> m_point_data;
    std::unordered_map<std::string, Field> m_cell_data;
};

template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, void>
VtrFile::point_data(const std::string& name, std::vector<std::vector<T>>& data) {
    m_point_data[name] = {.type=VtkType::get<T>(), .data={}};
    std::vector<std::byte>& d1 = m_point_data[name].data;
    d1.resize(m_x.size() * m_y.size() * sizeof(T));
    if (data.size() != m_x.size()) {
        throw std::runtime_error("Size mismatch 'x'");
    }
    int nx = m_x.size();
    T* out = reinterpret_cast<T*>(d1.data());
    for (int i = 0; i < m_x.size(); ++i) {
        if (data[i].size() != m_y.size()) {
            throw std::runtime_error("Size mismatch 'y'");
        }
        for (int j = 0; j < m_y.size(); ++j) {
            out[nx * j + i] = data[i][j];
        }
    }
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, void>
VtrFile::point_data(const std::string& name, std::vector<std::vector<std::vector<T>>>& data) {
    m_point_data[name] = {.type=VtkType::get<T>(), .data={}};
    std::vector<std::byte>& d1 = m_point_data[name].data;
    d1.resize(m_x.size() * m_y.size() * m_z.size() * sizeof(T));
    if (data.size() != m_x.size()) {
        throw std::runtime_error("Size mismatch 'x'");
    }
    int nx = m_x.size();
    int ny = m_y.size();
    T* out = reinterpret_cast<T*>(d1.data());
    for (int i = 0; i < m_x.size(); ++i) {
        if (data[i].size() != m_y.size()) {
            throw std::runtime_error("Size mismatch 'y'");
        }
        for (int j = 0; j < m_y.size(); ++j) {
            if (data[i][j].size() != m_z.size()) {
                throw std::runtime_error("Size mismatch 'z'");
            }
            for (int k = 0; k < m_z.size(); ++k) {
                out[nx * ny * k + nx * j + i] = data[i][j][k];
            }
        }
    }
}

} // namespace zephyr::io