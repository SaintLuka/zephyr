#pragma once

#include <memory>
#include <set>
#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/box.h>


namespace zephyr::mesh::decomp {

using zephyr::geom::Vector3d;
using zephyr::geom::Box;

/// @brief Простенькая взвешенная диаграмма Вороного.
/// Нет быстрых функций поиска разбиения, но есть пара функций
/// для графического отображения.
/// @details Написано плохо и неэффективно. Переписать надо,
/// или найти говорую библиотеку, может поручить студенту?
class VDiagram {
public:
    // Характеристики балансировки

    double mobility    = 0.2;   // Скорость движения генераторов
    double growth_rate = 0.02;  // Скорость изменения весов
    double centroidal  = 0.05;  // Влияние смещения к центру масс

    /// Фукнции инициализации

    /// @brief Конструктор по умолчанию
    VDiagram() = default;

    /// @brief Конструктор со случайными генераторами
    VDiagram(const Box& domain, int size);

    /// @brief Диаграмма по генераторам
    VDiagram(const Box& domain, const std::vector<Vector3d>& gs);


    /// Функции доступа

    /// @return Число ячеек диаграммы
    int size() const;

    /// @return Центр ячейки Воронного
    const Vector3d& coords(int idx) const;

    /// @return Координаты центра масс ячейки
    std::vector<Vector3d>& centers();

    /// @return Координаты x генераторов
    std::vector<double> coords_x() const;
    std::vector<double> coords_y() const;
    std::vector<double> coords_z() const;

    /// @return Координаты x центров масс ячеек
    std::vector<double> centers_x() const;
    std::vector<double> centers_y() const;
    std::vector<double> centers_z() const;

    /// @return Веса генераторов
    std::vector<double> weights() const;

    /// @return Степень вершин графа из генераторов ячеек Воронного и связей
    std::vector<int> degrees();

    /// @return Хроматическое число диаграммы
    /// @details На самом деле не хроматическое число, а число цветов,
    /// получившееся в алгоритме раскраски
    int chromatic_number() const;

    /// @return Цвета ячеек
    const std::vector<int>& colors();

    /// @return Список координат x границ ячеек
    std::vector<std::vector<double>> lines_x();

    /// @return Список координат y границ ячеек
    std::vector<std::vector<double>> lines_y();

    /// @return Соедининения смежных генераторов
    std::vector<std::vector<double>> connections_x();
    std::vector<std::vector<double>> connections_y();

    /// @return Список радиусов вписаных окружностей
    const std::vector<double>& search_radii();

    double search_radius(int iGen) const;

    /// @return Границы вписанных окружностей (координата x)
    std::vector<std::vector<double>> search_area_x();

    /// @return Границы вписанных окружностей (координата x)
    std::vector<std::vector<double>> search_area_y();

    /// @brief Установить координаты и вес генератора
    void add_generator(double x, double y, double w);

    /// @brief Установить генератор
    void set_coords(int iGen, double x, double y);

    /// @brief Установить вес генератора
    void set_weight(int iGen, double w);

    /// @brief Получить вес генератора
    double get_weight(int iGen);

    /// @brief Получить коодинату генератора
    double get_coord_x(int iGen);
    double get_coord_y(int iGen);
    double get_coord_z(int iGen);

    /// @brief Установить новые положения генераторов
    void set_coords(const std::vector<Vector3d>& coords);

    /// @brief Установить все веса диаграммы
    void set_weights(const std::vector<double>& ws);

    /// @brief Раскрасить диаграмму
    /// @param K Количество цветов
    void paint();


    /// Весовая диаграмма

    /// @brief Принадлежность ячейке диаграммы
    int rank(const Vector3d& v) const;

    /// @brief Расстояние от точки p до генератора g с весом w
    static double wdistance(const Vector3d& p, const Vector3d& g, double w = 0.0);

    /// @brief Расстояние от точки p до генератора iGen
    /// @param p точка
    /// @param iGen номер генератора
    double wdistance(const Vector3d& p, int iGen) const;

    /// @brief Расстояние от генератора i до генератора j
    /// @param i номер первого генератора
    /// @param j номер второго генератора
    double distance_gen(int i, int j);

    /// @brief Функция возвращает 0 на границе ячейки, положительное значение
    /// внутри и отрицательное снаружи
    /// @param p Точка
    /// @param iGen Номер генератора
    double edge_function(const Vector3d& p, int iGen) const;


    /// Функции балансировки

    /// @brief Балансировка по площадям
    void balancing();

    /// @brief Балансировка взвешенной диаграммы
    void balancing(const std::vector<double>& loads);

private:
    /// @brief Перевести диаграму в неактуальное состояние. Функция вызывается
    /// после смещения генераторов или изменения весов диаграммы.
    void changed();


    /// @brief Нормализовать веса диаграммы, сумма = 0.0
    void normalize();

    /// @brief Непосредственное построение диаграммы, поиск границ ячеек,
    /// поиск смежных ячеек
    void build();

    /// @brief Ограничивающий прямоугольник
    Box m_domain;

    /// @brief Координаты генераторов диаграммы Воронного
    std::vector<Vector3d> m_coords;

    /// @brief Веса генераторов
    std::vector<double> m_weights;



    /// @brief Истинно, если диаграмма полностью построена
    bool m_actual;

    /// @brief Координаты центров масс ячеек
    std::vector<Vector3d> m_centers;

    /// @brief Границы ячеек Воронного
    std::vector<std::vector<Vector3d>> m_lines;

    /// @brief Радиусы вписанных в ячейки окружностей с центрами в генераторах
    std::vector<double> m_search_radii;

    /// @brief Список смежных подобластей
    std::vector<std::set<int>> m_adjacency;

    /// @brief Цвета ячеек (для красивого отображения)
    std::vector<int> m_colors;
};

} // namespace zephyr::mesh::decomp