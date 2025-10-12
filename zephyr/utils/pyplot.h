#pragma once

#include <vector>
#include <string>
#include <variant>
#include <optional>
#include <iostream>

#include <zephyr/configuration.h>

#ifdef ZEPHYR_PYTHON
#define nopython
#else
#define nopython { std::cerr << "Enable ZEPHYR_PYTHON to use utils::pyplot functions\n"; }
#endif

namespace zephyr::utils {

/// @brief Опции отрисовки фигуры
struct figure_options {
    // figure size width x heigh in inches
    std::tuple<double, double> figsize = {-1.0, -1.0};
    std::optional<int>         dpi     = 130;
};

/// @brief Опции и стиль линий в plot, loglog и подобных
struct line_options {
    std::optional<std::string> linestyle  = std::nullopt;
    std::optional<double>      linewidth  = std::nullopt;
    std::optional<std::string> color      = std::nullopt;
    std::optional<std::string> marker     = std::nullopt;
    std::optional<double>      markersize = std::nullopt;
    std::optional<std::string> label      = std::nullopt;
};

/// @brief Опции заполнения полигона (pyplot.fill)
struct fill_options {
    std::optional<std::string> color= std::nullopt;
};

/// @brief Опции стрелочки pyplot.arrow
struct arrow_options {
    std::optional<std::string> face_color  = std::nullopt;
    std::optional<std::string> edge_color  = std::nullopt;
    std::optional<double>      head_length = std::nullopt;
    std::optional<double>      head_width  = std::nullopt;
};

/// @brief Опции фукнции scatter plot
struct scatter_options {
    /// @brief Размер маркера. Единственно число, либо массив чисел
    std::optional<std::variant<double, std::vector<double>>> s;

    /// @brief Цвета маркеров. Единственная строка-цвет или массив чисел.
    std::optional<std::variant<std::string, std::vector<double>>> c;

    std::optional<std::string> marker;  ///< Маркер
    std::optional<std::string> cmap;    ///< Используемый colormap

    /// @brief Отображение чисел в цветовую палитру, возможные варианты
    /// (их бы изучить подробнее):
    //      "linear": LinearScale,
    //      "log":    LogScale,
    //      "symlog": SymmetricalLogScale,
    //      "asinh":  AsinhScale,
    //      "logit":  LogitScale,
    //      "function": FuncScale,
    //      "functionlog": FuncScaleLog
    std::optional<std::string> norm;

    /// @brief Границы отображения
    std::optional<double> vmin, vmax;


    bool colorbar = false; ///< Отрисовать colorbar (для 3D)
    std::string clabel;    ///< Подпись для colorbar (для 3D)
};

/// @brief Опции отрисовки трёхмерной поверхности
struct surface_options {
    std::optional<std::string> cmap  = std::nullopt;
};

/// @brief Простая обёртка для matplotlib.pyplot
///
/// @details Использует непосредственно matplotlib.pyplot через pybind11.
/// Все зависимости спрятаны в .cpp файл. При первом вызове конструктора
/// создает python-интерпретатор, который висит до завершения программы.
class pyplot final {
public:
    using array_t   = std::vector<double>;
    using array2d_t = std::vector<array_t>;

    /// @brief При первом вызове конструктора создается python-интерпретатор,
    /// который существует до завершения программы.
    pyplot() nopython;

    ~pyplot() nopython;

    /// @{ Создание фигуры и осей

    /// @brief Создать Figure
    void figure(const figure_options& args = {}) const nopython;

    /// @brief Создать Figure с заданы м номером
    void figure(int num, const figure_options& args = {}) const nopython;

    /// @brief Классический subplot, оси нумеруются от 0 до nrows*ncols
    void subplot(int nrows, int ncols, int plot_idx) const nopython;

    /// @brief Улучшенный subplot, позволяет точнее выбрать место для осей
    void subplot2grid(int nrows, int ncols, int rowid, int colid) const nopython;

    void tight_layout() const nopython;

    void show() const nopython;

    /// @}

    /// @{ Настройки осей

    void grid(bool enable) const nopython;

    void set_aspect_equal() const;

    void xlim(double xmin, double xmax) const nopython;

    void ylim(double xmin, double xmax) const nopython;

    /// @}

    /// @{ Различные подписи

    void title(const std::string& title) const nopython;

    void suptitle(const std::string& title) const nopython;

    void legend() const nopython;

    void xlabel(const std::string& text) const nopython;

    void ylabel(const std::string& text) const nopython;

    void text(double x, double y, const std::string& text) const nopython;

    /// @}

    /// @{ Двумерные графики

    /// @brief Просто нарисовать точку
    void marker(double x, double y, const line_options& args = {}) const nopython;

    /// @brief Нарисовать стрелку
    void arrow(double x, double y, double dx, double dy, const arrow_options& args = {}) const nopython;

    /// @brief Классический plot
    void plot(const array_t& x, const array_t& y, const line_options& args = {}) const nopython;

    /// @brief Классический plot со строкой формата вроде format="b--".
    void plot(const array_t& x, const array_t& y, const char* format, const line_options& args = {}) const nopython;

    void semilogx(const array_t& x, const array_t& y, const line_options& args = {}) const nopython;
    void semilogx(const array_t& x, const array_t& y, const char* format, const line_options& args = {}) const nopython;

    void semilogy(const array_t& x, const array_t& y, const line_options& args = {}) const nopython;
    void semilogy(const array_t& x, const array_t& y, const char* format, const line_options& args = {}) const nopython;

    void loglog(const array_t& x, const array_t& y, const line_options& args = {}) const nopython;
    void loglog(const array_t& x, const array_t& y, const char* format, const line_options& args = {}) const nopython;

    /// @brief Нарисовать полигон
    void fill(const array_t& x, const array_t& y, const fill_options& args = {}) const nopython;

    /// @brief Построить множество точек
    void scatter(const array_t &x, const array_t &y, const scatter_options& args = {}) const nopython;

    /// @brief Находит какие-то оси и строит colorbar (для 3D осей не работает)
    void colorbar(const std::string& label = "") const nopython;

    /// @}

    /// @{ Трёхмерные графики

    /// @brief Создать 3D оси и трёхмерный scatter
    void scatter3D(const array_t &x, const array_t &y, const array_t &z, const scatter_options& args = {}) const nopython;

    /// @brief Создать 3D оси и построить поверхность
    void plot_surface(const array2d_t &x, const array2d_t &y, const array2d_t &z, const surface_options& args = {}) const nopython;

    /// @}

#ifdef ZEPHYR_PYTHON
protected:
    void plot_impl(const char* func_name, const array_t& x, const array_t& y, const line_options& args = {}) const nopython;

    void plot_impl(const char* func_name, const array_t& x, const array_t& y, const char* format, const line_options& args = {}) const nopython;

    class pyport;
    pyport* impl;
#endif
};

} // namespace zephyr::utils