#pragma once

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
/// или найти готовую библиотеку, может поручить студенту?
class VDiagram {
public:
    // Характеристики балансировки

    double mobility    = 0.2;   // Скорость движения генераторов
    double growth_rate = 0.02;  // Скорость изменения весов
    double centroidal  = 0.05;  // Влияние смещения к центру масс

    /// Функции инициализации

    /// @brief Конструктор по умолчанию
    VDiagram() = default;

    /// @brief Конструктор со случайными генераторами
    VDiagram(const Box& domain, int size);

    /// @brief Диаграмма по генераторам
    VDiagram(const Box& domain, const std::vector<Vector3d>& gs);


    /// Функции доступа

    /// @return Число ячеек диаграммы
    int size() const;

    /// @return Центр ячейки Вороного
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

    /// @return Степень вершин графа из генераторов ячеек Вороного и связей
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

    /// @return Соединения смежных генераторов
    std::vector<std::vector<double>> connections_x();
    std::vector<std::vector<double>> connections_y();

    /// @return Список радиусов вписанных окружностей
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

    /// @brief Установить генератор
    void set_coords(int iGen, const Vector3d& p);

    /// @brief Установить вес генератора
    void set_weight(int iGen, double w);

    /// @brief Получить вес генератора
    double get_weight(int iGen) const;

    /// @brief Получить координату генератора
    double get_coord_x(int iGen) const;
    double get_coord_y(int iGen) const;
    double get_coord_z(int iGen) const;

    Vector3d get_coord(int iGen) const;

    /// @brief Установить новые положения генераторов
    void set_coords(const std::vector<Vector3d>& coords);

    /// @brief Установить все веса диаграммы
    void set_weights(const std::vector<double>& ws);

    /// @brief Раскрасить диаграмму
    void paint();


    /// Весовая диаграмма

    /// @brief Принадлежность ячейке диаграммы
    int rank(const Vector3d& v) const;

    /// @brief Расстояние от точки p до генератора g с весом w
    static double wdistance(const Vector3d& p, const Vector3d& g, double w = 0.0);

    /// @brief Расстояние от точки p до генератора iGen
    /// @param p точка
    /// @param iGen Номер генератора
    double wdistance(const Vector3d& p, int iGen) const;

    /// @brief Расстояние от генератора i до генератора j
    /// @param i Номер первого генератора
    /// @param j Номер второго генератора
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
    /// @brief Перевести диаграмму в неактуальное состояние. Функция вызывается
    /// после смещения генераторов или изменения весов диаграммы.
    void changed();


    /// @brief Нормализовать веса диаграммы, сумма = 0.0
    void normalize();

    /// @brief Непосредственное построение диаграммы, поиск границ ячеек,
    /// поиск смежных ячеек
    void build();

    /// @brief Ограничивающий прямоугольник
    Box domain_;

    /// @brief Координаты генераторов диаграммы Вороного
    std::vector<Vector3d> coords_;

    /// @brief Веса генераторов
    std::vector<double> weights_;



    /// @brief Истинно, если диаграмма полностью построена
    bool actual_;

    /// @brief Координаты центров масс ячеек
    std::vector<Vector3d> centers_;

    /// @brief Границы ячеек Вороного
    std::vector<std::vector<Vector3d>> lines_;

    /// @brief Радиусы вписанных в ячейки окружностей с центрами в генераторах
    std::vector<double> search_radii_;

    /// @brief Список смежных подобластей
    std::vector<std::set<int>> adjacency_;

    /// @brief Цвета ячеек (для красивого отображения)
    std::vector<int> colors_;
};

} // namespace zephyr::mesh::decomp