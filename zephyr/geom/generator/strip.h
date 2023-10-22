#pragma once

#include <zephyr/configuration.h>
#include <zephyr/geom/primitives/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @class Strip. Простой класс для создания Storage для хранения
/// квазиодномерной сетки (прямоугольная сетка толщиной в одну ячейку).
class Strip : public Generator {
public:
    using Ptr = std::shared_ptr<Strip>;

    /// @brief Тип задания узлов
    enum class Type {
        UNIFORM,
        RANDOM
    };

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left  = Boundary::ZOE;
        Boundary right = Boundary::ZOE;
    };

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit Strip(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    Strip(double xmin, double xmax, Type nodes = Type::UNIFORM);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Strip::Ptr create(Args&&... args){
        return std::make_shared<Strip>(std::forward<Args>(args)...);
    }

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    void set_nx(int nx);

    /// @brief Установить желаемое число ячеек сетки
    void set_size(int N);

    /// @brief Установить флаги граничных условий
    void set_boundaries(Boundaries bounds);

    /// @brief Количество ячеек сетки
    int size() const final;

    /// @brief Ограничивающий объем
    Box bbox() const final;

    /// @brief Создать сетку общего вида
    Grid make() final;


    // Далее не самые полезные get-функции

    /// @brief Левая граница
    double x_min() const;

    /// @brief Правая граница
    double x_max() const;

    /// @brief Нижняя граница
    double y_min() const;

    /// @brief Правая граница
    double y_max() const;

    /// @brief Число ячеек по оси X
    int nx() const;

    /// @brief Есть ли периодичность по оси X?
    bool periodic_along_x() const;

private:
    /// @brief Проверить параметры сетки перед созданием
    void check_params() const final;

    ///@brief Соотношение сторон прямоугольника
    const double aspect = 0.1;

    /// @brief Тип задания узлов сетки
    Type m_type;

    /// @brief Число ячеек сетки
    int m_nx;

    /// @brief Левая и правая граница полосы
    double m_xmin, m_xmax;

    /// @brief Граничные условия слева и справа
    Boundaries m_bounds;
};

} // namespace zephyr::geom::generator
