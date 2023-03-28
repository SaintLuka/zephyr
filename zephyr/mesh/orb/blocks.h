#pragma once

#include <array>
#include <memory>

#include <zephyr/data/type/Vector3d.h>
#include <zephyr/network/decomposition/ORB/transform.h>

namespace zephyr { namespace network { namespace decomposition {

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

    /// @brief Конструктор с полным или частичным заданием размеров декомозиции.
    /// Для одномерного разбиения достаточно указать n_proc.
    /// Для двумерного разбиения требуется указать n_proc и n_proc_1.
    /// Для трехмерного разбиения требуется указать n_proc, n_proc_1 и n_proc_2.
    /// @param type Тип декомпозиции
    /// @param n_proc Полное число блоков
    /// @param n_proc_1, n_proc_2, n_proc_3 Число блоков вдоль осей, порядок
    /// соответствует порядку переменных в типе декомпозиции. В качестве
    /// неопределенных значений следует использовать отрицательные числа.
    Blocks(const std::string& type, size_t n_proc,
            int n_proc_1, int n_proc_2, int n_proc_3);

    // Автоматический конструктор по геометрическим размерам.
    // Внимание: не работает с полярным разбиением.
    Blocks(const std::string& type, size_t n_proc, const Vector3d& sizes);

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
    size_t size() const;

    /// @brief Число блоков в колонке
    size_t size(int i) const;

    /// @brief Число блоков в колонке и строке
    size_t size(int i, int j) const;

    /// @brief Центр блока
    Vector3d center(int rank) const;

    /// @brief Находится ли точка v с процесса my_rank в R-окрестности
    /// процесса neib_rank
    bool is_near(int my_rank, int neib_rank, const Vector3d& v, double R);

    /// @brief Балансировка нагрузки
    void balancing_simple(const std::vector<double>& loads, double alpha);

    /// @brief Балансировка нагрузки
    void balancing_newton(const std::vector<double>& loads, double alpha);

    /// @brief Границы блоков
    std::vector<std::vector<Vector3d>> lines() const;

    /// @brief Вывести информацию о структуре
    void info() const;

private:
    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Одномерная декомпозиция
    void init_sizes_1D();

    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Двумерная декомпозиция
    void init_sizes_2D(int n_proc_1);

    /// @brief Проинициализировать массивы m_nx, m_ny, m_nz
    /// @details Трехмерная декомпозиция
    void init_sizes_3D(int n_proc_1, int n_proc_2);

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

    /// @brief Границы кубоида
    Vector3d m_vmin, m_vmax;

    /// @brief Число блоков
    size_t m_size;
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

} // decomposition
} // network
} // zephyr