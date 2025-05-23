#pragma once

#include <functional>

#include <zephyr/mesh/euler/amr_storage.h>

namespace zephyr::mesh {

class Children {
public:
    Children(AmrStorage::Iterator ch1, AmrStorage::Iterator ch2,
             AmrStorage::Iterator ch3, AmrStorage::Iterator ch4)
            : data({ch1, ch2, ch3, ch4,
                    nullptr, nullptr, nullptr, nullptr}) {}

    Children(AmrStorage::Iterator ch1, AmrStorage::Iterator ch2,
             AmrStorage::Iterator ch3, AmrStorage::Iterator ch4,
             AmrStorage::Iterator ch5, AmrStorage::Iterator ch6,
             AmrStorage::Iterator ch7, AmrStorage::Iterator ch8)
            : data({ch1, ch2, ch3, ch4,
                    ch5, ch6, ch7, ch8}) {}

    Children(const std::array<AmrStorage::Iterator, 4>& arr)
            : data({arr[0], arr[1], arr[2], arr[3],
                    nullptr, nullptr, nullptr, nullptr}) {}

    Children(std::array<AmrStorage::Iterator, 8>&& arr)
            : data(std::move(arr)) {}

    int count() const {
        return data[4] ? 8 : 4;
    }

    int datasize() const {
        return data[0].datasize();
    }

    AmrStorage::Item &operator[](int idx) {
        return *data[idx];
    };

    struct Iterator {

        Iterator(AmrStorage::Iterator *ptr)
                : m_ptr(ptr) {}

        AmrStorage::Item &operator*() {
            return *(*m_ptr);
        }

        void operator++() {
            ++m_ptr;
        }

        bool operator!=(const Iterator &it) const {
            return m_ptr != it.m_ptr;
        }

    private:
        AmrStorage::Iterator *m_ptr;
    };

    Iterator begin() {
        return {data.data()};
    }

    Iterator end() {
        return {data.data() + count()};
    }


private:
    std::array<AmrStorage::Iterator, 8> data;
};

class QCell;
class SoaChildren;

/// @brief Класс содержит набор функций для огрубления и распределения
// (физических) данных в ячейках при адаптации. Определяются пользователем.
struct Distributor {
    /// @brief Тип функции, распределяющей данные между дочерними ячейками
    using split_function = std::function<void(AmrStorage::Item&, Children &)>;

    /// @brief Тип функции, объединяющей данные дочерних ячеек
    using merge_function = std::function<void(Children &, AmrStorage::Item&)>;

    /// @brief Тип функции, распределяющей данные между дочерними ячейками
    using split_function_soa = std::function<void(QCell&, SoaChildren&)>;

    /// @brief Тип функции, объединяющей данные дочерних ячеек
    using merge_function_soa = std::function<void(SoaChildren&, QCell&)>;


    split_function split;  ///< Распределение данным между дочерними
    merge_function merge;  ///< Объединение данных дочерних ячеек

    split_function_soa split_soa;  ///< Распределение данным между дочерними
    merge_function_soa merge_soa;  ///< Объединение данных дочерних ячеек

    /// @brief Конструктор по умолчанию. Определяет Distributor, который
    /// вообще ничего не делает
    Distributor();

    /// @brief Создает Distributor, который вообще ничего не делает.
    static Distributor empty();

    /// @brief Создает Distributor, который определяет функции split2D, split3D,
    /// merge2D, merge3D простейшим способом. Для функций split используется
    /// простой перенос данных в дочерние ячейки. Для функций merge
    /// используется перенос данных в родительсую ячейку из первой дочерней.
    static Distributor simple();

};

} // namespace zephyr::mesh