#pragma once

#include <filesystem>
#include <zephyr/utils/range.h>

#include <zephyr/geom/side.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/euler/amr_verts.h>
#include <zephyr/mesh/euler/amr_faces.h>

// forward declaration для классов из geom
namespace zephyr::geom {
struct Box;
class Quad;
class SqQuad;
class Cube;
class SqCube;
class Line;
class Polygon;
class Polyhedron;
class Grid;

namespace generator {
class Strip;
class Rectangle;
class Cuboid;
}
}

namespace zephyr::mesh {

/// @brief Опции сетки
struct MeshOpts {
    int  dim      = -1;    ///< Размерность сетки
    bool adaptive = true;  ///< Адаптивная сетка?
    bool linear   = true;  ///< Линейная адаптивная сетка?
    bool axial    = false; ///< Осевая симметрия
    bool nodes    = false; ///< Уникальные узлы
};

/// @brief Набор ячеек в форме Structure of Arrays (набор массивов).
///
/// Поддерживается три типа сеток:
///   1. Двумерная AMR сетка, по 8 граней, по 9 вершин на ячейку.
///   2. Трехмерная AMR сетка, по 24 грани, по 27 вершин на ячейку.
///   3. Неструктурированная/произвольная сетка. Произвольное число граней
///      и вершин на ячейку, но вершины не уникальны.
///
/// Сетка не поддерживает уникальные вершины и грани. Таким образом, одна
/// внутренняя грань сетки хранится в двух экземплярах. Члены класса (за
/// исключением базовых) помечены как public. Доступ к ним открыт, как
/// если бы это была обычная структура. Сеточные данные обрабатываются
/// специальными методами.
class AmrCells final {
    // aliases inside class
    using Vector3d = geom::Vector3d;

    /// @brief Характеристическая функция (функция-индикатор)
    using InFunction = std::function<bool(const Vector3d &)>;

    /// @brief Пространственная функция
    using SpFunction = std::function<double(const Vector3d &)>;

    /// @{ @name Глобальные характеристики ячеек

    index_t m_size  = 0;      ///< Число ячеек

    int  dim_      = -1;     ///< Размерность ячейки
    bool adaptive_ = false;  ///< Адаптивная ячейка?
    bool linear_   = true;   ///< Линейная ячейка?
    bool axial_    = false;  ///< Осевая симметрия?

    /// @}

public:
    /// @{ @name Основные индексы ячейки

    std::vector<int>     rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)
    std::vector<index_t> index;  ///< Глобальный индекс элемента в локальном Storage
                                 /// (< 0 для неопределенных элементов, элементов на удаление)
    /// @}
    /// @{ @name Характеристики адаптивных ячеек

    std::vector<int> flag;       ///< Желаемый флаг адаптации
    std::vector<int> level;      ///< Уровень адаптации (0 для базовой)
    std::vector<index_t> b_idx;  ///< Индекс среди базовых ячеек
    std::vector<index_t> z_idx;  ///< Индекс ячейки на z-кривой

    /// @}
    /// @{ @name Геометрия ячеек

    std::vector<Vector3d> center;      ///< Барицентр ячейки
    std::vector<double>   volume;      ///< Объем трехмерной / площадь двумерной ячейки
    std::vector<double>   volume_alt;  ///< Объем двумерной осесимметричной ячейки

    /// @}
    /// @{ @name Грани и вершины ячеек

    AmrFaces faces;  ///< Массив граней ячеек
    AmrVerts verts;  ///< Массив вершин ячеек

    /// @}

    Storage  data;   ///< Данные ячеек

public:
    /// @{ @name Конструкторы

    /// @brief Базовый конструктор
    /// @param options Настройки сетки
    explicit AmrCells(MeshOpts options = {});

    /// @brief Пустое множество ячеек с таким же набором опций и типов
    AmrCells same() const;

    /// @}

    /// @{ @name Общие характеристики ячеек

    /// @brief Размерность сетки
    int dim() const { return dim_; }

    /// @brief Сетка допускает адаптацию?
    bool adaptive() const { return adaptive_; }

    /// @brief Сетка с осевой симметрией?
    bool axial() const { return axial_; }

    /// @brief Используются линейные AMR-ячейки (или квадратичные)
    bool linear() const { return linear_; }

    /// @brief Сетка хранит уникальные узлы?
    bool has_nodes() const { return verts.has_nodes(); }

    /// @brief Настройки сетки
    MeshOpts options() const;

    /// @}

    /// @{ @name Размеры хранилища

    /// @brief Пустое хранилище?
    bool empty() const { return m_size == 0; }

    /// @brief Число ячеек (синоним)
    index_t n_cells() const { return m_size; }

    /// @brief Полное число граней
    index_t n_faces() const { return faces.n_faces(); }

    /// @brief Полное число вершин с дубликатами
    index_t n_verts() const { return verts.n_verts(); }

    /// @brief Очистить хранилище
    void clear();

    /// @brief Только для AMR ячеек с установленной размерностью
    void resize_amr(index_t n_cells);

    /// @brief Только для AMR ячеек с установленной размерностью
    void reserve_amr(index_t n_cells);

    /// @brief Увеличить массивы под ячейки, грани и вершины
    /// @param n_cells Полное число ячеек
    /// @param n_faces Полное число граней
    /// @param n_nodes Полное число вершин
    void resize(index_t n_cells, index_t n_faces, index_t n_nodes);

    /// @brief Зарезервировать под неструктурированные ячейки
    /// @param n_cells Полное число ячеек
    /// @param n_faces Полное число граней
    /// @param n_nodes Полное число вершин
    void reserve(index_t n_cells, index_t n_faces, index_t n_nodes);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Расход памяти (только массивы ячеек)
    memory_t memory_usage() const;

    /// @}

    /// @{ @name Топологические свойства ячеек

    /// @brief Актуальная ячейка?
    bool is_actual(index_t ic) const { return index[ic] >= 0; }

    /// @brief Ячейка к удалению
    bool is_undefined(index_t ic) const { return index[ic] < 0; }

    /// @brief Устанавливает index = -1 (ячейка вне сетки)
    void set_undefined(index_t ic) { index[ic] = -1; }

    /// @brief Число актуальных граней ячейки, для адаптивной ячейки может
    /// быть меньше faces.max_count, для неструктурированной ячейки (полигон
    /// или многогранник совпадает с faces.max_count)
    int face_count(index_t ic) const;

    /// @brief Название грани AMR-ячейки
    std::string face_name(index_t ic, index_t iface) const {
        return geom::side_to_string(iface - faces.offsets[ic], dim_);
    }

    /// @}

    /// @{ @name Геометрические свойства ячеек

    /// @brief Линейный размер ячейки по оси x (от левой до правой грани)
    double hx(index_t ic) const;

    /// @brief Линейный размер ячейки по оси y (от нижней до верхней грани)
    double hy(index_t ic) const;

    /// @brief Линейный размер ячейки по оси z (от задней до передней грани)
    double hz(index_t ic) const;

    /// @brief Линейный размер ячейки
    double linear_size(index_t ic) const {
        return dim_ < 3 ? std::sqrt(volume[ic]) : std::cbrt(volume[ic]);
    }

    /// @brief Обычный объем или объем осесимметичной ячейки
    double get_volume(index_t ic, bool axial) const {
        return axial ? volume_alt[ic] : volume[ic];
    }

    /// @brief Диаметр вписанной окружности.
    /// @details Для AMR-ячейки представляет собой минимальное расстояние между
    /// противоположными гранями. Для полигона --- диаметр вписанной окружности
    /// для правильного многоугольника аналогичной площади.
    /// Величину удобно использовать совместно с условием Куранта.
    /// Для двумерных расчетов на прямоугольных сетках совпадает с минимальной
    /// стороной прямоугольной ячейки.
    double incircle_diameter(index_t ic) const;

    /// @brief Получить вершину по индексу внутри ячейки
    const Vector3d& vertex(index_t ic, int iv) const {
        return verts[verts.offsets[ic] + iv];
    }

    /// @brief Bounding box ячейки
    geom::Box bbox(index_t ic) const;

    /// @brief Создать полигон из ячейки (для 3D ячеек -- UB)
    geom::Polygon polygon(index_t ic) const;

    /// @brief Создать полигон из ячейки (для 2D ячеек -- UB)
    geom::Polyhedron polyhedron(index_t ic) const;

    /// @}

    /// @{ @name Сечения и интегрирование по ячейке

    /// @brief Оценка объемной доли, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @details Относительно быстрая функция, проверяет функцию inside только
    /// в узлах ячейки, позволяет быстро выяснить, содержит ли ячейка
    /// границу двух областей. Если ячейка внутри тела, то возвращает строго
    /// единицу 1.0, если снаружи -- строго ноль 0.0.
    double approx_vol_fraction(index_t ic, const InFunction& inside) const;

    /// @brief Объемная доля, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @param n_points Число тестовых точек, для которых проверяется функция
    /// inside, погрешность определения объемной доли ~ 1/N.
    double volume_fraction(index_t ic, const InFunction& inside, int n_points) const;

    /// @brief Функция func является константой в ячейке?
    /// @details Проверяется значение функции в узлах и в центре ячейки,
    /// если все значения совпадают, то считается, что функция принимает
    /// постоянное значение в пределах ячейки.
    bool const_function(index_t ic, const SpFunction& func) const;

    /// @brief Интеграл скалярной функции по ячейке
    /// @param n_points Разбиение по сторонам
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(index_t ic, const SpFunction& func, int n_points) const;

    /// @}

    /// @{ @name Работа с данными

    void move_item(index_t from, index_t to);

    /// @brief Скопировать все данные с индекса from на индекс to.
    void copy_data(index_t from, index_t to);

    /// @brief Скопировать все данные целиком с индекса from,
    /// в хранилище dst на индекс to
    void copy_data(index_t from, AmrCells* dst, index_t to) const;

    /// @brief Скопировать ячейку с позиции ic в хранилище cells на индекс jc,
    /// грани на позицию iface, вершины на позицию inode.
    void copy_geom(index_t ic, AmrCells& cells,
            index_t jc, index_t face_beg, index_t node_beg) const;

    void copy_geom_basic(index_t ic, AmrCells& cells,
            index_t jc, index_t face_beg, index_t node_beg) const;

    /// @}

    /// @{ @name Конструирование и добавление ячеек

    /// @brief Двумерная простая
    void set_cell(index_t ic, const geom::Quad &quad);

    /// @brief Двумерная с осевой симметрией
    void set_cell(index_t ic, const geom::Quad &quad, bool axial);

    /// @brief Двумерная криволинейная
    void set_cell(index_t ic, const geom::SqQuad &quad);

    /// @brief Двумерная криволинейная с осевой симметрией
    void set_cell(index_t ic, const geom::SqQuad &quad, bool axial);

    /// @brief Трехмерная простая
    void set_cell(index_t ic, const geom::Cube &cube);

    /// @brief Трехмерная криволинейная ячейка
    void set_cell(index_t ic, const geom::SqCube &cube);

    /// @brief Добавить простую двумерную ячейку в виде отрезка.
    /// Сохраняется как четырехугольник
    void push_back(const geom::Line &line);

    /// @brief Добавить в конец двумерную полигональную ячейку.
    /// Не адаптивная ячейка, может представлять четырехугольник, но при этом
    /// вершины будут упорядочены иначе.
    void push_back(const geom::Polygon &poly);

    /// @brief Добавить в конец трехмерную ячейку-многогранник.
    /// Не адаптивная ячейка, может представлять шестигранник в виде куба,
    /// но при этом вершины будут упорядочены иначе.
    void push_back(const geom::Polyhedron& poly);

    /// @}

    /// @{ @name Функции для дебага

    /// @brief Вывести информацию о ячейке
    void print_info(index_t ic) const;

    /// @brief Вывести информацию о ячейке в виде python скрипта
    /// для визуализации
    void visualize(index_t ic, std::string filename) const;

    /// @brief Проверить базовую геометрию ячейки
    /// @return -1 для плохой ячейки
    int check_geometry(index_t ic) const;

    /// @brief Проверить ориентацию граней
    /// @return -1 для плохой ячейки
    int check_base_face_orientation(index_t ic) const;

    /// @brief Проверить порядок вершин
    /// @return -1 для плохой ячейки
    int check_base_vertices_order(index_t ic) const;

    /// @brief Проверить сложные грани
    /// @return -1 для плохой ячейки
    int check_complex_faces(index_t ic) const;

    /// @brief Проверка связности ячеек для однопроцессорной версии
    int check_connectivity(index_t ic) const;

    /// @brief Проверка связности ячеек в MPI версии
    int check_connectivity(index_t ic, const AmrCells& ghosts) const;

    /// @brief Полное сохранение сетки
    /// @param root Корневая директория для бэкапа (существует и пустая)
    /// @param file Открытый файл для записи (backup.json)
    /// @param tab Отступ секции в json файле
    /// @param variables Список имен переменных для сохранения
    void backup(const std::filesystem::path& root, std::ofstream& file,
        const std::string& tab, const std::vector<std::string>& variables) const;

    /// @}

protected:
    /// @brief Увеличить только массивы данных ячеек
    void resize_cells(index_t n_cells);

    /// @brief Увеличить буферы только для массивов данных ячеек
    void reserve_cells(index_t n_cells);

    /// @brief Увеличить буферы только для массивов данных ячеек
    void shrink_to_fit_cells();

    void push_back_impl(const geom::Polyhedron& poly);
};

} // namespace zephyr::mesh