#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/maps.h>

namespace zephyr::geom {

/// @brief Базовый класс для описания двумерного полигона
/// При наличии массива вершин можно создать полигон на его основе, при
/// этом будет использован тот же буффер памяти без дополнительного
/// копирования. К примеру,
/// std::vector<Vector3d> verts;    // Создаем массив и заполняем
///                                    некоторым образом
/// Polygon(verts).sort();          // Упорядочивает массив вершин
///                                    против часовой стрелки
///                                    (без выделения памяти)
/// auto S = Polygon(verts).area(); // Вычисляем площадь полигона
///                                    (без выделения памяти)
///
/// Другой вариант: использование дочерних классов PolygonS или PolygonD
/// (лучше данный класс, если операция создания полигона является массовой
/// или вызывается на тредах)
class Polygon {
public:

    /// @brief Создание полигона на основе имеющегося массива вершин
    /// заданого размера
    /// @warning Допускается изменение массива, запрещено удалять массив
    /// до уничножения экземпляра Polygon.
    Polygon(std::vector<Vector3d>& vertices);

    /// @brief Копирование запрещено
    Polygon(const Polygon& ) = delete;

    /// @brief Перемещение запрещено
    Polygon(Polygon&& ) = delete;

    /// @brief Число вершин в буфере
    inline int size() const { return m_size; };

    /// @brief Пустой полигон?
    inline bool empty() const { return m_size < 3; }

    /// @brief Оператор доступа по индексу
    inline Vector3d& operator[](int idx) { return vs[idx]; }

    /// @brief Оператор доступа по индексу
    inline const Vector3d& operator[](int idx) const { return vs[idx]; }

    /// @brief Центр полигона (среднее вершин)
    Vector3d center() const;

    /// @brief Сортировать вершины по часовой стрелке
    /// @param c Точка внутри полигона, если известена
    void sort(const Vector3d& c = {NAN, NAN, NAN});

    /// @brief Точка внутри полигона?
    bool inside(const Vector3d& p) const;

    /// @brief Площадь произвольного многоугольника
    /// Вершины должны располагаться в плоскости (x, y).
    /// @param c Точка внутри полигона, если известена
    double area(const Vector3d& c = {NAN, NAN, NAN}) const;

    /// @brief Барицентр произвольного многоугольника.
    /// Вершины должны располагаться в плоскости (x, y).
    /// @param area Площадь многоугольника (если известна)
    Vector3d centroid(double area = 0.0) const;

    /// @brief Объем фигуры вращения (вокруг оси x).
    /// Вершины должны располагаться в плоскости (x, y).
    double volume_as() const;

protected:
    /// @brief Создание полигона на основе имеющегося массива вершин
    /// заданого размера
    /// Используется в конструкторах PolygonS и PolygonD
    Polygon(Vector3d* buff, int size);

    struct MinMax {
        Vector3d min;
        Vector3d max;
    };

    /// @brief Минимальная и максимальныя (по проекции на ось n)
    /// вершины многоугольника
    MinMax minmax(const Vector3d& n) const;

    /// @brief Отсечь от полигона часть прямой с внешней нормалью n,
    /// проходящей через точку p
    /// @param part Выходной параметр, многоугольник, который отсекается
    void clip(const Vector3d& p, const Vector3d& n, Polygon& part) const;

    /// @brief Отсечь от полигона часть прямой с внешней нормалью n,
    /// проходящей через точку p
    /// @param part Выходной параметр, многоугольник, который отсекается
    /// @param slice Выходной параметр, сечение.
    void clip(const Vector3d& p, const Vector3d& n, Polygon& part, Line& slice) const;

    /// @brief Площадь многоугольника, отсекаемая прямой с внешней
    /// нормалью n, проходящей через точку p
    /// @param part Выходной параметр, многоугольник, который отсекается
    double clip_area(const Vector3d& p, const Vector3d& n, Polygon& part) const;

public:
    /// @brief Площадь и длина сечения
    struct AnS { double area, slice; };

    /// @brief Площадь многоугольника, отсекаемая прямой с внешней
    /// нормалью n, проходящей через точку p
    /// @param part Выходной параметр, многоугольник, который отсекается
    AnS clip_area_and_slice(const Vector3d& p, const Vector3d& n, Polygon& part) const;

protected:
    /// @brief Находит отсечение от полигона с заданой объемной долей
    /// @param part Выходной параметр, многоугольник, который отсекается
    Vector3d find_section(const Vector3d& n, double alpha, Polygon& part) const;

    /// @brief Метод дихотомии
    Vector3d find_section_bisect(const Vector3d& n, double alpha, Polygon& part) const;

    /// @brief Метод Ньютона
    Vector3d find_section_newton(const Vector3d& n, double alpha, Polygon& part) const;


    int m_size;    ///< Число вершин полигона
    Vector3d* vs;  ///< Указатель на первую вершину массива
};

template <int max_size = 8>
class PolygonS : public Polygon {
public:
    /// @brief Контруктор по умолчанию, создает пустой полигон
    PolygonS() : Polygon(buffer.data(), 0) { }

    /// @brief Конструктор копирования
    PolygonS(const PolygonS<max_size>& poly)
        : Polygon(buffer.data(), poly.size()),
          buffer(poly.buffer) {
    }

    /// @brief Простых путей мы не ищем, вместо инициализации с использованием
    /// initializer_list делаем рекурсивное заполнение по переменному числу
    /// входных параметров. Создание класса:
    /// PolygonS poly = {a, b, c, ... }, где a, b, c - имеют тип Vector3d
    /// Число аргументов должно быть равно max_size.
    template<class...Args, typename = typename
            std::enable_if<sizeof...(Args) == max_size>::type>
    PolygonS(Args... args) : Polygon(buffer.data(), sizeof...(Args)) {
        fill(buffer.data(), std::forward<Args>(args)...);
    }

    /// @brief Изменить размер полигона
    void resize(int new_size) {
        m_size = std::min(new_size, max_size);
    }

    /// @brief Создание четырехугольника из уголовых точек
    /// адаптивной ячейки.
    PolygonS(const Quad& quad) : Polygon(buffer.data(), 4) {
        vs[0] = quad.vs<-1, -1>();
        vs[1] = quad.vs<+1, -1>();
        vs[2] = quad.vs<+1, +1>();
        vs[3] = quad.vs<-1, +1>();
    }

    /// @brief Отсечь от полигона часть прямой с внешней нормалью n,
    /// проходящей через точку p
    PolygonS<max_size + 1> clip(const Vector3d& p, const Vector3d& n) const {
        PolygonS<max_size + 1> part;

        Polygon::clip(p, n, part);
        return part;
    }

    /// @brief Площадь многоугольника, отсекаемая прямой с внешней
    /// нормалью n, проходящей через точку p
    double clip_area(const Vector3d& p, const Vector3d& n) const {
        PolygonS<max_size + 1> part;
        return Polygon::clip_area(p, n, part);
    }

    /// @brief Находит отсечение от полигона с заданой объемной долей
    /// @param part Выходной параметр, многоугольник, который отсекается
    Vector3d find_section(const Vector3d& n, double alpha) const {
        PolygonS<max_size + 1> part;
        return Polygon::find_section(n, alpha, part);
    }

protected:

    /// @brief Заполнение элементов массива по одному
    template <class... Args>
    void fill(Vector3d* ptr, const Vector3d& x, Args... args) {
        *ptr = x;
        fill(ptr + 1, std::forward<Args>(args)...);
    }

    /// @brief Заполнение элементов массива по одному (база рекурсии)
    void fill(Vector3d* ptr, const Vector3d& x) {
        *ptr = x;
    }

    /// @brief Массив вершин
    std::array<Vector3d, max_size> buffer;
};

using PolyTri  = PolygonS<3>;
using PolyQuad = PolygonS<4>;

class PolygonD : public Polygon {
public:
    PolygonD(const std::vector<Vector3d>& verts)
        : Polygon(nullptr, 0) {
        buffer = verts;
        update();
    }

    PolygonD(int size)
            : Polygon(nullptr, size), buffer(size) {
        update();
    }

    void reserve(int size) {
        buffer.reserve(size);
        update();
    }

    void resize(int size) {
        buffer.resize(size);
        update();
    }

private:
    void update() {
        vs = buffer.data();
        m_size = buffer.size();
    }

    std::vector<Vector3d> buffer;
};

} // namespace zephyr::geom