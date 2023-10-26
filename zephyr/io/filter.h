#pragma once

#include <memory>
#include <set>

#include <zephyr/mesh/storage.h>

namespace zephyr { namespace io {

using zephyr::mesh::AmrStorage;

/// @class Абстрактный фильтр для сохранения элементов хранилища в VTU файл
class Filter {
public:
    /// @brief Оператор фильтрации, возвращает true, если элемент хранилища,
    /// переданный в качетсве аргумента, необходимо сохранить.
    virtual bool operator()(AmrStorage::Item&) const = 0;

    /// @brief Является ли оператор фильтрации тривиальным.
    /// Тривиальным будем считать оператор фильтрации, который оставляет все
    /// элементы хранилища.
    virtual bool is_trivial() const = 0;
};

/// @class Тривиальный фильтр, сохраняет все элементы хранилища
class TrivialFilter : public Filter {
public:
    bool operator()(AmrStorage::Item& cell) const override {
        return true;
    }

    bool is_trivial() const override {
        return true;
    }
};

/// @class Оставляет элементы только внутри некоторого параллелепипеда
class BoxFilter : public Filter {
public:
    BoxFilter(
        double xmin, double xmax,
        double ymin, double ymax,
        double zmin = -std::numeric_limits<double>::infinity(),
        double zmax = +std::numeric_limits<double>::infinity())
        : xmin(xmin), xmax(xmax),
          ymin(ymin), ymax(ymax),
          zmin(zmin), zmax(zmax) {
    }

    bool operator()(AmrStorage::Item& cell) const override {
        auto& c = cell.center;
        return xmin <= c.x() && c.x() <= xmax &&
               ymin <= c.y() && c.y() <= ymax &&
               zmin <= c.z() && c.z() <= zmax;
    }

    bool is_trivial() const override {
        return false;
    }

private:
    double xmin, xmax, ymin, ymax, zmin, zmax;
};

/// @class Оставляет только элементы внутри некоторого цилиндра
class CylinderFilter : public Filter {
public:
    CylinderFilter(double r, double x = 0.0, double y = 0.0)
        : r(r), x(x), y(y) {
    }

    bool operator()(AmrStorage::Item& cell) const override {
        auto& c = cell.center;
        return (c.x() - x) * (c.x() - x) + (c.y() - y) * (c.y() - y) <= r * r;
    }

    bool is_trivial() const override {
        return false;
    }

private:
    double r, x, y;
};

/// @class Сложный фильтр, позволяет объединить набор других фильтров.
/// При пустом наборе фильтров является тривиальным. Новые фильтры добавляются
/// с помощью изящного оператора +=.
/// @code
/// ComplexFilter filter;
/// filter += CellsFilter();
/// filter += BoxFilter(0.0, 1.0, 0.0, 1.0);
/// @endcode
class ComplexFilter : public Filter {
public:
    bool operator()(AmrStorage::Item& cell) const override {
        for (auto& filter: filters) {
            if (!(*filter)(cell)) {
                return false;
            }
        }
        return true;
    }

    bool is_trivial() const override {
        return filters.empty();
    }

    /// @brief Добавить новый фильтр к составному
    template<class SomeFilter>
    ComplexFilter& append(const SomeFilter& filter) {
        if (!filter.is_trivial()) {
            static_assert(std::is_base_of<Filter, SomeFilter>::value, "Add strange filter");
            filters.emplace_back(new SomeFilter(filter));
        }
        return *this;
    }

    /// @brief Добавить новый фильтр к составному
    template<class SomeFilter>
    ComplexFilter& operator+=(const SomeFilter& filter) {
        if (!filter.is_trivial()) {
            static_assert(std::is_base_of<Filter, SomeFilter>::value, "Add strange filter");
            filters.emplace_back(new SomeFilter(filter));
        }
        return *this;
    }

    /// @brief Сбросить все фильтры
    void reset() {
        filters.clear();
    }

private:
    std::vector<std::unique_ptr<Filter>> filters;
};

} // namespace io
} // namespace zephyr

