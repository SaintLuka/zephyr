#pragma once

#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/cell_type.h>
#include <zephyr/mesh/index.h>

namespace zephyr::mesh {

/// @brief Индексы смежности граней
/// @details Короткое объяснение. Если сосед через грань
///             с этого процесса  |  с другого процесса
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    alien :    < 0             |    < aliens.size
struct AmrAdjacent final {

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    std::vector<int> rank;

    /// @brief Индекс смежной ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    std::vector<index_t> index;

    /// @brief Индекс смежной ячейки в массиве aliens (или -1, если соседняя
    /// ячейка с данного процесса).
    std::vector<index_t> alien;

    /// @brief Индекс базовой ячейки, которая содержит грань
    std::vector<index_t> basic;

    /// @brief Расширить массивы по числу граней
    void resize(index_t n_faces);

    /// @brief Расширить массивы по числу граней
    void reserve(index_t n_faces);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Локальная соседняя ячейка?
    bool is_local(index_t iface) const { return alien[iface] < 0; }

    /// @brief Удаленная соседняя сосед?
    bool is_alien(index_t iface) const { return alien[iface] >= 0; }

    /// @brief Получить хранилище ячеек, в котором находится сосед, а также
    /// индекс соседа в данном хранилище
    template<class SomeArray>
    std::tuple<const SomeArray &, index_t> get_neib(index_t iface,
            const SomeArray &locals, const SomeArray &aliens) const {
        if (alien[iface] < 0) {
            return {locals, index[iface]};
        } else {
            return {aliens, alien[iface]};
        }
    }
};

/// @brief Перечисление используется для выбора граней
/// с определенным направлением нормалей
enum class Direction : int {
    ANY = 0,  // Любое направление нормали
    X   = 1,
    Y   = 2,
    Z   = 3
};

/// @brief Грани AMR-сетки, развернутые в SoA.
struct AmrFaces final {
    // aliases inside class
    using Boundary = zephyr::geom::Boundary;
    using Vector3d = zephyr::geom::Vector3d;

    AmrAdjacent adjacent; ///< Индексы смежных ячеек

    std::vector<Boundary> boundary;  ///< Тип граничного условия
    std::vector<Vector3d> normal;    ///< Внешняя нормаль к грани
    std::vector<Vector3d> center;    ///< Барицентр грани
    std::vector<double>   area;      ///< Площадь грани
    std::vector<double>   area_alt;  ///< "Альтернативная" площадь грани

    /// @brief Список индексов вершин в массиве вершин ячейки
    std::vector<std::array<int, 4>> vertices;


    /// @brief Число граней
    index_t size() const { return boundary.size(); }

    /// @brief Изменить размер под число граней
    void resize(index_t n_faces);

    /// @brief Расширить буффер под число граней
    void reserve(index_t n_faces);

    /// @brief Сжать до актуальных размеров
    void shrink_to_fit();

    /// @brief Добавить грани, соответствующие ячейке
    /// @details Есть реализация вставки для разных типов ячеек.
    void insert(index_t iface, geom::CellType ctype, int count = -1);

    /// @brief Является ли грань граничной?
    bool is_boundary(index_t iface) const;

    /// @brief Является ли грань актуальной?
    bool is_actual(index_t iface) const;

    /// @return 'true', если грань не актуальна
    bool is_undefined(index_t iface) const;

    /// @brief Установить неопределенную грань
    void set_undefined(index_t iface);

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает с выбраным направлением
    bool to_skip(index_t iface, Direction dir) const;

    /// @brief Внешняя нормаль грани на площадь
    Vector3d area_n(index_t iface) const;

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double get_area(index_t iface, bool axial = false) const;

    /// @brief Симметричная точка относительно грани
    Vector3d symm_point(index_t iface, const Vector3d& p) const;
};


inline bool AmrFaces::is_boundary(index_t iface) const {
    return boundary[iface] != Boundary::ORDINARY &&
           boundary[iface] != Boundary::PERIODIC &&
           boundary[iface] != Boundary::UNDEFINED;
}

inline bool AmrFaces::is_actual(index_t iface) const {
    return boundary[iface] != Boundary::UNDEFINED;
}

inline bool AmrFaces::is_undefined(index_t iface) const {
    return boundary[iface] == Boundary::UNDEFINED;
}

inline void AmrFaces::set_undefined(index_t iface) {
    boundary[iface] = Boundary::UNDEFINED;
    adjacent.rank[iface]  = -1;
    adjacent.index[iface] = -1;
    adjacent.alien[iface] = -1;
    adjacent.basic[iface] = -1;
}

inline AmrFaces::Vector3d AmrFaces::area_n(index_t iface) const {
    return area[iface] * normal[iface];
}

inline double AmrFaces::get_area(index_t iface, bool axial) const {
    return axial ? area_alt[iface] : area[iface];
}

inline AmrFaces::Vector3d AmrFaces::symm_point(index_t iface, const Vector3d& p) const {
    return p + 2.0 * (center[iface] - p).dot(normal[iface]) * normal[iface];
}

} // namespace zephyr::mesh