#pragma once

#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Простой класс для создания AmrStorage для хранения
/// квазиодномерной сетки (прямоугольная сетка шириной в одну ячейку).
class Strip final : public Generator {
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

    /// @brief Конструктор класса по конфигурации
    explicit Strip(const Json& config);

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    /// @param nodes Тип распределения узлов
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

    /// @brief Осевая симметрия странно
    void set_axial(bool) override { m_axial = false; }

    /// @brief Нелинейная сетка странно
    void set_linear(bool) override { m_linear = true; }

    /// @brief Ограничивающий объем
    Box bbox() const override;

    /// @brief Создать сетку общего вида
    Grid make() const override;

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

    /// @brief Граничные условия
    Boundaries bounds() const { return m_bounds; }

    /// @brief Есть ли периодичность по оси X?
    bool periodic_along_x() const;

private:
    /// @brief Проверить параметры сетки перед созданием
    void check_params() const;

    ///@brief Соотношение сторон прямоугольника
    const double aspect = 0.01;

    /// @brief Тип задания узлов сетки
    Type m_type;

    /// @brief Число ячеек сетки
    int m_nx{0};

    /// @brief Левая и правая граница полосы
    double m_xmin, m_xmax;

    /// @brief Граничные условия слева и справа
    Boundaries m_bounds;
};

} // namespace zephyr::geom::generator
