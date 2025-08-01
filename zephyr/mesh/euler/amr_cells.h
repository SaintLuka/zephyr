#pragma once

#include <zephyr/utils/range.h>
#include <zephyr/utils/storage.h>
#include <zephyr/mesh/primitives/side.h>
#include <zephyr/mesh/euler/amr_faces.h>

// Forward Declaration для классов из geom
namespace zephyr::geom {
class Quad;
class SqQuad;
class Cube;
class SqCube;
class Line;
class Polygon;
class Polyhedron;
}

namespace zephyr::mesh {

/// @brief Квадратичное отображение на квадрат/куб в зависимости от размерности
template <int dim>
using SqMap = std::conditional_t<dim < 3, geom::SqQuad, geom::SqCube>;

/// @brief Аналог AmrCell развернутый в структуру массивов,
/// теоретически, все функции тоже можно просто скопировать.
/// Три типа сеток.
/// 1. Двумерная AMR сетка, по 8 граней, по 9 вершин на ячейку.
/// 2. Трехмерная AMR сетка, по 24 грани, по 27 вершин на ячейку.
/// 3. Неструктурированная/произвольная сетка. Произвольное число
/// граней и вершин на ячейку, но вершины не уникальны.
/// Нет неактуальных вершин.
/// Неактуальные грани допускаются для AMR ячеек.
///
/// Почему всё public, обоснование
class AmrCells final {
    // aliases inside class
    using Vector3d = geom::Vector3d;
    using SoaStorage = utils::SoaStorage;
    using AmrVerts = std::vector<Vector3d>;

    index_t m_size  = 0;      ///< Число ячеек

    int  m_dim      = -1;     ///< Размерность ячейки
    bool m_adaptive = false;  ///< Адаптивная ячейка?
    bool m_linear   = true;   ///< Линейная ячейка?
    bool m_axial    = false;  ///< Осевая симметрия?

public:

    // Главные индексы ячейки

    std::vector<int>     rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)
    std::vector<index_t> index;  ///< Глобальный индекс элемента в локальном Storage
                                 /// ( < 0 для неопределенных элементов, элементов на удаление)

    // Величины, связанные с адаптацией

    std::vector<int> flag;       ///< Желаемый флаг адаптации
    std::vector<int> level;      ///< Уровень адаптации (0 для базовой)
    std::vector<index_t> b_idx;  ///< Индекс среди базовых ячеек
    std::vector<index_t> z_idx;  ///< Индекс ячейки на z-кривой

    // Геометрия ячеек

    std::vector<Vector3d> center;      ///< Барицентр ячейки
    std::vector<double>   volume;      ///< Объем трехмерной / площадь двумерной ячейки
    std::vector<double>   volume_alt;  ///< Объем двумерной осесимметричной ячейки

    // Грани и вершины ячеек

    std::vector<index_t> face_begin;   ///< Индекс первой грани ячейки
    std::vector<index_t> node_begin;   ///< Индекс первой вершины ячейки

    AmrFaces   faces;  ///< Массив граней ячеек
    AmrVerts   verts;  ///< Массив вершин ячеек

    SoaStorage data;   ///< SoA с данными ячеек



    // Конструкторы

    /// @brief Никак не уйти, иногда нужны
    AmrCells() = default;

    /// @brief Базовый конструктор
    /// @param dim Размерность сетки (2 для осевой симметрии)
    /// @param adaptive Использовать возможность адаптации?
    /// @param axial Сетка с осевой симметрией?
    explicit AmrCells(int dim, bool adaptive = false, bool axial = false);


    /// @brief Создать пустое хранилище с таким же набором типов
    AmrCells same() const;


    /// @brief Размерность сетки
    int dim() const { return m_dim; }

    /// @brief Сетка допускает адаптацию?
    bool adaptive() const { return m_adaptive; }

    /// @brief Сетка с осевой симметрией?
    bool axial() const { return m_axial; }

    /// @brief Используются линейные AMR-ячейки (или квадратичные)
    bool linear() const { return m_linear; }

    /// @brief Изменить размерность
    void set_dimension(int dim);

    /// @brief Использовать адаптивные ячейки
    void set_adaptive(bool adaptive = true);

    /// @brief Использовать осевую симметрию
    void set_axial(bool axial = true);

    /// @brief Использовать линейные отображения
    void set_linear(bool linear);



    // Базовые свойства всего хранилища

    /// @brief Пустое хранилище?
    bool empty() const { return m_size == 0; }

    /// @brief Число ячеек (не очевидно, что речь о ячейках)
    index_t size() const { return m_size; }

    /// @brief Число ячеек (синоним)
    index_t n_cells() const { return m_size; }

    /// @brief Полное число граней
    index_t n_faces() const { return faces.size(); }

    /// @brief Полное число вершин
    index_t n_nodes() const { return verts.size(); }


    /// @{ @name Топологические свойства ячеек

    /// @brief Актуальная ячейка?
    bool is_actual(index_t ic) const { return index[ic] >= 0; }

    /// @brief Ячейка к удалению
    bool is_undefined(index_t ic) const { return index[ic] < 0; }

    /// @brief Устанавливает index = -1 (ячейка вне сетки)
    void set_undefined(index_t ic) { index[ic] = -1; }

    /// @brief Число актуальных граней ячейки, для адаптивной ячейки может
    /// быть меньше max_face_count, для неструктурированной ячейки (полигон
    /// или многогранник совпадает с max_face_count)
    int face_count(index_t ic) const;

    /// @brief Максимальное число граней для ячейки
    int max_face_count(index_t ic) const {
        return face_begin[ic + 1] - face_begin[ic];
    }

    /// @brief Полный диапазон граней ячейки (могут встречаться неактуальные)
    utils::range<index_t> faces_range(index_t ic) const {
        return utils::range(face_begin[ic], face_begin[ic + 1]);
    }

    /// @brief Число вершин, оно же максимальное, хранение неактуальных вершин
    /// сейчас не допускается.
    int node_count(index_t ic) const {
        return node_begin[ic + 1] - node_begin[ic];
    }

    /// @brief Число вершин, оно же максимальное, хранение неактуальных вершин
    /// сейчас не допускается.
    int max_node_count(index_t ic) const {
        return node_begin[ic + 1] - node_begin[ic];
    }

    /// @brief Полный диапазон вершин ячейки
    utils::range<index_t> nodes_range(index_t ic) const {
        return utils::range(node_begin[ic], node_begin[ic + 1]);
    }

    /// @brief Простая грань на стороне?
    template <int dim>
    bool simple_face(index_t ic, Side<dim> side) const {
        return faces.is_undefined(face_begin[ic] + side[1]);
    }

    /// @brief Сложная грань на стороне?
    template <int dim>
    bool complex_face(index_t ic, Side<dim> side) const {
        return faces.is_actual(face_begin[ic] + side[1]);
    }

    /// @brief Название грани для AMR-ячейки
    std::string face_name(index_t ic, index_t iface) const {
        return side_to_string(iface - face_begin[ic], m_dim);
    }

    /// @}


    /// @{ @name Геометрические свойства ячеек

    /// @brief Линейный размер ячейки
    double linear_size(index_t ic) const {
        return m_dim < 3 ? std::sqrt(volume[ic]) : std::cbrt(volume[ic]);
    }

    /// @brief Обычный объем или объем осесимметичной ячейки
    double get_volume(index_t ic, bool axial) const {
        return axial ? volume_alt[ic] : volume[ic];
    }

    /// @brief Диаметр вписаной окружности.
    /// @details Для AMR-ячейки представляет собой минимальное расстояние между
    /// противоположными гранями. Для полигона --- диаметр вписаной окружности
    /// для правильного многоугольника аналогичной площади.
    /// Величину удобно использовать совместно с условием Куранта.
    /// Для двумерных расчетов на прямоугольных сетках совпадает с минимальной
    /// стороной прямоугольной ячейки.
    double incircle_diameter(index_t ic) const;

    /// @brief Указатель на первую вершину
    Vector3d* vertices_data(index_t ic) {
        return verts.data() + node_begin[ic];
    }

    /// @brief Константный указатель на первую вершину
    const Vector3d* vertices_data(index_t ic) const {
        return verts.data() + node_begin[ic];
    }

    /// @brief Ссылка на вешины в форме набора узлов квадратичного отображения
    template <int dim>
    SqMap<dim>& mapping(index_t ic) {
        return *reinterpret_cast<SqMap<dim>*>(vertices_data(ic));
    }

    /// @brief Ссылка на вешины в форме набора узлов квадратичного отображения
    template <int dim>
    const SqMap<dim>& mapping(index_t ic) const {
        return *reinterpret_cast<const SqMap<dim>*>(vertices_data(ic));
    }

    /// @brief Создать полигон из ячейки (для 3D ячеек -- UB)
    geom::Polygon polygon(index_t ic) const;


    /// @{ @name Сечения и интегрирование по ячейке

    /// @brief Оценка объемной доли, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @details Относительно быстрая функция, проверяет функцию inside только
    /// на узлах ячейки, позволяет быстро выяснить, содержит ли ячейка
    /// границу двух областей. Если ячейка внутри тела, то возвращает строго
    /// единицу 1.0, если снаружи -- строго ноль 0.0.
    double approx_vol_fraction(index_t ic, const std::function<double(const Vector3d &)> &inside) const;

    /// @brief Объемная доля, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @param n_points Число тестовых точек, для которых проверяется функция
    /// inside, погрешность определения объемной доли ~ 1/N.
    double volume_fraction(index_t ic, const std::function<double(const Vector3d &)> &inside, int n_points) const;

    /// @brief Функция func является константой на ячейке?
    /// @details Проверяется значение функции в узлах и в центре ячейки,
    /// если все значения совпадают, то считается, что функция принимает
    /// постоянное значение в пределах ячейки.
    bool const_function(index_t ic, const std::function<double(const Vector3d&)>& func) const;

    /// @brief Интеграл скалярной функции по ячейке
    /// @param n_points Разбиение по сторонам
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(index_t ic, const std::function<double(const Vector3d&)>& func, int n_points) const;

    /// @}


    // Работа с данными

    void move_item(index_t ic);

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


    // Размеры хранилища

    /// @brief Очистить хранилище
    void clear();

    /// @brief Только для AMR ячеек с установленной размерностью
    void resize_amr(index_t n_cells);

    /// @brief Только для AMR ячеек с установленной размерностью
    void reserve_amr(index_t n_cells);

    /// @brief Увеличить массивы под ячейки, грани и вершины
    /// @param n_faces Полное число граней!
    /// @param n_nodes Полное число вершин!
    void resize(index_t n_cells, index_t n_faces, index_t n_nodes);

    /// @brief Зарезервировать под неструктурированные ячейки
    /// @param n_faces Максимальное число граней на ячейку
    /// @param n_nodes Максимальное число вершин на ячейку
    void reserve(index_t n_cells, index_t n_faces, index_t n_nodes);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();



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
    int check_connectivity(index_t ic, const AmrCells& aliens) const;

    /// @}


protected:
    /// @brief Увеличить только массивы данных ячеек
    void resize_cells(index_t n_cells);

    /// @brief Увеличить буфферы только для массивов данных ячеек
    void reserve_cells(index_t n_cells);

    /// @brief Увеличить буфферы только для массивов данных ячеек
    void shrink_to_fit_cells();
};

} // namespace zephyr::mesh