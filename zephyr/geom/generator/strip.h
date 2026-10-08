#pragma once

#include <zephyr/configuration.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Простой класс для создания квазиодномерной сетки
/// (прямоугольная сетка шириной в одну ячейку).
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
    /// @param x_min, x_max Границы прямоугольника по оси x
    /// @param nodes Тип генерации узлов
    Strip(double x_min, double x_max, Type nodes = Type::UNIFORM);

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
    void set_axial(bool) override { axial_ = false; }

    /// @brief Нелинейная сетка странно
    void set_linear(bool) override { linear_ = true; }

    /// @brief Ограничивающий объем
    Box bbox() const override;

    /// @brief Создать сетку общего вида
    Grid make() const override;

    /// @brief Может инициализировать хранилище
    bool can_make_cells() const override { return true; }

    /// @brief Инициализация SoA-хранилища сетки
    mesh::RawCells make_cells(bool unique_nodes) const override;

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
    Boundaries bounds() const { return bounds_; }

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
    int nx_{0};

    /// @brief Левая и правая граница полосы
    double x_min_, x_max_;

    /// @brief Граничные условия слева и справа
    Boundaries bounds_;
};

} // namespace zephyr::geom::generator
