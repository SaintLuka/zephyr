#pragma once

#include <array>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/mesh/decomp/orb/transform.h>

namespace zephyr::mesh::decomp {

using zephyr::geom::Box;

struct barrier {
    double val;
    double dx;

    barrier() : val(0.0), dx(0.0) { }

    barrier(double _val) : val(_val), dx(0.0) { }

    operator double& () { return val; }

    operator const double& () const { return val; }
};

class Blocks {
public:
    /// @brief Тривиальный конструктор.
    /// Создает единственный блок на всё пространство.
    Blocks();

    /// @brief Конструктор с автоматической декомпозицией по размерам области.
    /// Подбирается так, чтобы блоки были близки к квадратам/кубам.
    /// @param domain Вычислительная область
    /// @param type Тип декомпозиции
    /// @param size Общее число блоков ( > 0 )
    Blocks(const Box& domain, const std::string& type, int size);

    /// @brief Конструктор с автоматической декомпозицией по размерам области.
    /// Подбирается так, чтобы блоки были близки к квадратам/кубам.
    /// @param domain Вычислительная область
    /// @param type Тип декомпозиции
    /// @param size Общее число блоков. Для одномерной декомпозиции можно
    /// указать size = -1, в этом случае параметр nx определяет размер
    /// @param nx Число блоков по первой координате
    Blocks(const Box& domain, const std::string& type, int size, int nx);

    /// @brief Конструктор с автоматической декомпозицией по размерам области.
    /// Подбирается так, чтобы блоки были близки к квадратам/кубам.
    /// @param domain Вычислительная область
    /// @param type Тип декомпозиции
    /// @param size Общее число блоков. Для двумерной декомпозиции можно
    /// указать size = -1, в этом случае параметр ny определяет размер
    /// @param nx Число блоков по первой координате
    Blocks(const Box& domain, const std::string& type, int size,
            const std::vector<int>& ny);

    /// @brief Прямое отображение
    Vector3d mapping(const Vector3d& ) const;

    /// @brief Обратное отображение
    Vector3d inverse(const Vector3d& ) const;

    /// @brief Установить границы
    void setup_bounds(Vector3d vmin, Vector3d vmax);

    /// @brief Обновить границы
    void update_bounds(const Vector3d& vmin, const Vector3d& vmax);

    /// @brief Ранг процесса
    /// @param v Вектор в локальных координатах
    int rank(const Vector3d& v) const;

    /// @brief Общее число блоков
    int size() const;

    /// @brief Число блоков в колонке
    int size(int i) const;

    /// @brief Число блоков в колонке и строке
    int size(int i, int j) const;

    /// @brief Центр блока
    Vector3d center(int rank) const;

    /// @brief Балансировка нагрузки
    void balancing_simple(const std::vector<double>& loads, double alpha);

    /// @brief Балансировка нагрузки
    void balancing_newton(const std::vector<double>& loads, double alpha);

    /// @brief Границы блоков
    std::vector<std::vector<Vector3d>> lines() const;

    /// @brief Вывести информацию о структуре
    void info() const;

private:
    /// @brief Инициализировать единичный блок
    void init_single();

    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Одномерная декомпозиция
    void init_sizes_1D(int nx);

    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Двумерная декомпозиция
    void init_sizes_2D(const std::vector<int>& ny);

    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Трехмерная декомпозиция
    void init_sizes_3D(const std::vector<std::vector<int>>& nz);

    /// @brief Проинициализировать пределы
    void init_limits();

    /// @brief Проставить ранги и индексы
    void init_indices();

    /// @brief vmin, vmax, а также размеры установлены,
    /// устанавливает начальные координаты
    void init_coords();

    /// @brief Координаты установлены, вычислить ширину
    void init_widths();

    /// @brief Инициализировать массивы с нагрузкой
    void init_loads();

    /// @brief Установить нагрузку
    void set_loads(const std::vector<double>& loads);

    /// @brief
    void readiness() const;

    static void move_simple(double alpha, const std::vector<double>& loads,
            const std::vector<double>& widths, std::vector<barrier>& coords);

    static void move_newton(double alpha, const std::vector<double>& loads,
            const std::vector<double>& widths, std::vector<barrier>& coords);


    /// @brief Готовность
    bool m_ready;

    /// @brief Отображение в кубоид
    Transform m_map;

    /// @brief Кубоид
    Box m_box;

    /// @brief Число блоков
    int m_size;
    int m_nx;
    std::vector<int> m_ny;
    std::vector<std::vector<int>> m_nz;

    /// @brief Координаты блоков по рангу
    std::vector<int> m_i;
    std::vector<int> m_j;
    std::vector<int> m_k;

    /// @brief Ранги блоков
    std::vector<std::vector<std::vector<int>>> m_ranks;

    /// @brief Границы блоков
    std::vector<barrier> m_x_coords;
    std::vector<std::vector<barrier>> m_y_coords;
    std::vector<std::vector<std::vector<barrier>>> m_z_coords;

    /// @brief Ширина блоков
    std::vector<double> m_x_widths;
    std::vector<std::vector<double>> m_y_widths;
    std::vector<std::vector<std::vector<double>>> m_z_widths;

    /// @brief Нагрузка
    std::vector<double> m_x_loads;
    std::vector<std::vector<double>> m_y_loads;
    std::vector<std::vector<std::vector<double>>> m_z_loads;
};

} // namespace zephyr::mesh::decomp