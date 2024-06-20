#pragma once

#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/bface.h>

#include <zephyr/mesh/euler/amr_storage.h>
#include <zephyr/mesh/euler/eu_face.h>


namespace zephyr::mesh {

using zephyr::geom::Adjacent;
using zephyr::geom::Vector3d;
using zephyr::geom::Vector2d;

/// @brief Ячейка (итератор для доступа к элементам класса EuMesh)
/// @details Копирует функции AmrStorage::Item, а также добавляет
/// специфические функции для ячейки: volume, center, faces.
class EuCell {
public:
    // Конструкторы

    /// @brief Конструктор копирования
    EuCell(const EuCell &cell) = default;

    /// @brief Конструктор перемещения
    EuCell(EuCell &&cell) = default;

    /// @brief Основной конструктор. Взять ячейку из хранилища по индексу.
    /// Ячейка может быть из locals или из aliens. Если индекс меньше
    /// locals.size, тогда ячейка берется из хранилища locals, иначе
    /// выбирается ячейка aliens[idx - locals.size]. Если индекс слишком
    /// большой (больше locals.size + aliens.size), тогда возвращается
    /// итератор на конец хранилища locals.
    /// @param idx Индекс ячейки.
    EuCell(AmrStorage &locals, AmrStorage &aliens, int idx);

    /// @brief Вспомогательный конструктор. Взять смежную ячейку из locals
    /// или aliens хранилища по индексу смежности с грани.
    /// @param adj Индекс смежности ячейки.
    EuCell(AmrStorage &locals, AmrStorage &aliens, const Adjacent &adj);


    // Функции итератора

    /// @brief Разыменование итератора (для цикла for)
    inline EuCell &operator*() {
        return *this;
    }

    /// @brief Следующая ячейка
    inline EuCell &operator++() {
        ++m_it;
        return *this;
    }

    /// @brief Ячейка через step
    inline EuCell &operator+=(int step) {
        m_it += step;
        return *this;
    }

    /// @brief Ячейка через step
    inline EuCell operator+(int step) {
        EuCell res(*this);
        res += step;
        return res;
    }

    /// @brief Расстояние между двумя ячейками
    int operator-(const EuCell& cell) const {
        return m_it - cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator<(const EuCell& cell) {
        return m_it < cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator!=(const EuCell& cell) {
        return m_it != cell.m_it;
    }


    // Получение данных, обход граней и соседей

    /// @brief Указатель на данные ячейки
    inline Byte* data() {
        return m_it->data();
    }

    /// @brief Указатель на данные ячейки
    inline const Byte* data() const {
        return m_it->data();
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    inline U& operator()(const U &) {
        return *reinterpret_cast<U *>(m_it->data());
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    inline const U& operator()(const U &) const {
        return *reinterpret_cast<const U *>(m_it->data());
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    inline U& data(const U &) {
        return *reinterpret_cast<U *>(m_it->data());
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    inline const U& data(const U &) const {
        return *reinterpret_cast<const U *>(m_it->data());
    }

    /// @brief Ссылка на конкретные данные ячейки
    template<class T>
    inline T& operator()(const VarExtra<T> & var) {
        return *reinterpret_cast<T *>(m_it->data() + var.offset);
    }

    /// @brief Ссылка на конкретные данные ячейки
    template<class T>
    inline const T& operator()(const VarExtra<T> & var) const {
        return *reinterpret_cast<const T *>(m_it->data() + var.offset);
    }

    /// @brief Ссылка на геометрию ячейки
    inline geom::AmrCell& geom(){
        return *m_it;
    }

    /// @brief Ссылка на геометрию ячейки
    inline const geom::AmrCell& geom() const {
        return *m_it;
    }

    /// @brief Выбрать грань на стороне. Корректно работает
    /// для декартовых сеток без адаптации
    EuFace face(const geom::Side side) const;

    /// @brief Итератор по граням
    /// @param dir Выбрать грани по некоторым направлениям
    EuFaces faces(Direction dir = Direction::ANY) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(const BFace &face) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(const EuFace &face) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(int face_idx) const;

    /// @brief Получить указатель на данные соседа, минуя всякие промежуточные
    /// действия (должно работать быстро). Проводит проверки грани, если грань
    /// не определена или содержит граничное условие -- возвращается указатель
    /// на данные самой ячейки.
    const Byte* neib_data(const BFace& face) const;

    /// @brief Для данной ячейки вернуть ячейку из local AmrStorage
    /// @param idx Индекс ячейки в локальном хранилище
    EuCell locals(int idx);


    // Далее простые get/set функции

    /// @brief Размерность ячейки
    inline const int& dim() const { return m_it->dim; }

    /// @brief Индекс среди базовых ячеек
    inline const int& b_idx() const { return m_it->b_idx; }

    /// @brief Индекс ячейки на z-кривой
    inline const int& z_idx() const { return m_it->z_idx; }

    /// @brief Индекс новой ячейки (в алгоритмах)
    inline const int& next() const { return m_it->next; }

    /// @brief Уровень адаптации ячейки (0 для базовой)
    inline const int& level() const { return m_it->level; }

    /// @brief Желаемый флаг адаптации
    inline const int& flag() const { return m_it->flag; }

    /// @brief Желаемый флаг адаптации
    inline void set_flag(int flag) { m_it->flag = flag; }

    /// @brief Барицентр ячейки
    inline const Vector3d& center() const { return m_it->center; }

    /// @brief Барицентр ячейки
    inline const double& x() const { return m_it->center.x(); }

    /// @brief Барицентр ячейки
    inline const double& y() const { return m_it->center.y(); }

    /// @brief Барицентр ячейки
    inline const double& z() const { return m_it->center.z(); }

    /// @brief Линейный размер ячейки
    inline const double& size() const { return m_it->size; }

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    inline double volume() const { return m_it->volume(); }

    /// @brief Вершина ячейки по индексу
    inline const Vector3d& vs(int idx) const {
        return m_it->vertices[idx];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i, j, k in {-1, 0, +1}
    template <int i, int j, int k = -1>
    inline const Vector3d &vs() const {
        return m_it->vertices.vs<i, j, k>();
    }

    /// @brief Скорпировать вершины в полигон (двумерные ячейки)
    /// Для нелинейных AMR-ячеек возвращает до 8 граней.
    inline geom::Polygon polygon() const {
        return m_it->polygon();
    }

    /// @brief Диаметр вписаной окружности.
    /// @details Для AMR-ячейки представляет собой минимальное расстояние между
    /// противоположными гранями. Для полигона --- диаметр вписаной окружности
    /// для правильного многоугольника аналогичной площади.
    /// Величину удобно использовать совместно с условием Куранта.
    /// Для двумерных расчетов на прямоугольных сетках совпадает с минимальной
    /// стороной прямоугольной ячейки.
    inline double incircle_radius() const {
        return m_it->incircle_radius();
    }

    /// @brief Оценка объемной доли, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @details Относительно быстрая функция, проверяет функцию inside только
    /// на узлах ячейки, позволяет быстро выяснить, содержит ли ячейка
    /// границу двух областей. Если ячейка внутри тела, то возвращает строго
    /// единицу 1.0, если снаружи -- строго ноль 0.0.
    inline double approx_vol_fraction(const std::function<double(const Vector3d &)> &inside) const {
        return m_it->approx_vol_fraction(inside);
    }

    /// @brief Объемная доля, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @param n_points Число тестовых точек, для которых проверяется функция
    /// inside, погрешность определения объемной доли ~ 1/N.
    inline double volume_fraction(const std::function<double(const Vector3d &)> &inside, int n_points) const {
        return m_it->volume_fraction(inside, n_points);
    }


    // Функции для дебага

    /// @brief Вывести полную информаци о ячейке и её соседях
    void print_neibs_info() const;

    /// @brief Вывести полную информацию о ячейке
    inline void print_info() const { m_it->print_info(); }

    /// @brief Вывести информацию о ячейке в виде python скрипта
    inline void visualize() const { m_it->print_info(); }

protected:
    friend class EuFaces;

    AmrStorage::Iterator m_it;  ///< Итератор по хранилищу
    AmrStorage &m_locals;       ///< Ссылка на хранилище текущего процесса
    AmrStorage &m_aliens;       ///< Ссылка на хранилище соседей
};

} // namespace zephyr::mesh