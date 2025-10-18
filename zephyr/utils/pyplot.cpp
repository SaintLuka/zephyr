#include <iostream>

#include <zephyr/utils/pyplot.h>

#ifdef ZEPHYR_PYTHON
#include <pybind11/embed.h>
#include <pybind11/numpy.h>

namespace zephyr::utils {

namespace py = pybind11;

namespace {

// Двумерный массив в двумерный numpy array.
py::array_t<double> to_numpy(const pyplot::array2d_t& vec2d) {
    if (vec2d.empty()) {
        return py::array_t<double>(std::array<size_t, 2>{0, 0});
    }

    size_t rows = vec2d.size();
    size_t cols = vec2d[0].size();

    // Проверка размеров
    for (size_t i = 1; i < rows; ++i) {
        if (vec2d[i].size() != cols) {
            throw std::runtime_error("Все строки должны быть одинаковой длины");
        }
    }

    // Создаем плоский вектор для копирования
    std::vector<double> flat_data;
    flat_data.reserve(rows * cols);

    for (const auto& row : vec2d) {
        flat_data.insert(flat_data.end(), row.begin(), row.end());
    }

    // Создаем numpy массив из буфера
    return py::array_t<double>(
        {rows, cols},           // форма
        {cols * sizeof(double), sizeof(double)}, // шаги
        flat_data.data()        // данные
    );
}

py::kwargs get_kwargs(const figure_options& opts) {
    py::kwargs kwargs;
    if (std::get<0>(opts.figsize) >= 0.0 || std::get<1>(opts.figsize) >= 0.0) {
        kwargs["figsize"] = py::make_tuple(std::get<0>(opts.figsize), std::get<1>(opts.figsize));
    }
    if (opts.dpi.has_value()) { kwargs["dpi"] = py::int_(opts.dpi.value()); }
    return kwargs;
}

py::kwargs get_kwargs(const line_options& opts) {
    py::kwargs kwargs;
    if (opts.linestyle.has_value()) { kwargs["linestyle"] = py::str(opts.linestyle.value()); }
    if (opts.linewidth.has_value()) { kwargs["linewidth"] = py::float_(opts.linewidth.value()); }
    if (opts.color.has_value()) { kwargs["color"] = py::str(opts.color.value()); }
    if (opts.marker.has_value()) { kwargs["marker"] = py::str(opts.marker.value()); }
    if (opts.markersize.has_value()) { kwargs["markersize"] = py::float_(opts.markersize.value()); }
    if (opts.label.has_value()) { kwargs["label"] = py::str(opts.label.value()); }
    return kwargs;
}

py::kwargs get_kwargs(const fill_options& opts) {
    py::kwargs kwargs;
    if (opts.color.has_value()) { kwargs["color"] = py::str(opts.color.value()); }
    return kwargs;
}

py::kwargs get_kwargs(const arrow_options& opts) {
    py::kwargs kwargs;
    if (opts.face_color.has_value()) { kwargs["fc"] = py::str(opts.face_color.value()); }
    if (opts.edge_color.has_value()) { kwargs["ec"] = py::str(opts.edge_color.value()); }
    if (opts.head_length.has_value()) { kwargs["head_length"] = py::float_(opts.head_length.value()); }
    if (opts.head_width.has_value()) { kwargs["head_width"] = py::float_(opts.head_width.value()); }
    return kwargs;
}

py::kwargs get_kwargs(const scatter_options& args) {
    py::kwargs kwargs;
    if (args.s.has_value()) {
        if (std::holds_alternative<double>(args.s.value())) {
            kwargs["s"] = py::float_(std::get<double>(args.s.value()));
        }
        else {
            auto arr = std::get<std::vector<double>>(args.s.value());
            kwargs["s"] = py::array_t(arr.size(), arr.data());
        }
    }
    if (args.c.has_value()) {
        if (std::holds_alternative<std::string>(args.c.value())) {
            kwargs["c"] = py::str(std::get<std::string>(args.c.value()));
        }
        else {
            auto arr = std::get<std::vector<double>>(args.c.value());
            kwargs["c"] = py::array_t(arr.size(), arr.data());
        }
    }
    if (args.marker.has_value()) { kwargs["marker"] = py::str(args.marker.value()); }
    if (args.cmap.has_value()) { kwargs["cmap"] = py::str(args.cmap.value()); }
    if (args.norm.has_value()) { kwargs["norm"] = py::str(args.norm.value()); }
    if (args.vmin.has_value()) { kwargs["vmin"] = py::float_(args.vmin.value()); }
    if (args.vmax.has_value()) { kwargs["vmax"] = py::float_(args.vmax.value()); }
    return kwargs;
}

py::kwargs get_kwargs(const surface_options& opts) {
    py::kwargs kwargs;
    if (opts.cmap.has_value()) { kwargs["cmap"] = py::str(opts.cmap.value()); }
    return kwargs;
}

}

class pyplot::pyport {
public:
    py::module plt;
    py::module mpl_toolkits;
    py::module pylab_helpers;

    pyport() {
        try {
            if (!Py_IsInitialized()) {
                py::initialize_interpreter();
            }
            plt = py::module::import("matplotlib.pyplot");
            mpl_toolkits  = py::module::import("mpl_toolkits.mplot3d");
            pylab_helpers = py::module::import("matplotlib._pylab_helpers");

        } catch (const std::exception& e) {
            std::cerr << "Failed to initialize Python: " << e.what() << "\n";
        }
    }

    /// @brief Добавить трёхмерные оси на Figure
    py::object axes3D() const {
        py::object fig = plt.attr("gcf")(); // Текущая фигура
        return fig.attr("add_subplot")(111, py::arg("projection") = "3d");
    }
};

pyplot::pyplot() {
    impl = new pyport();
}

pyplot::~pyplot() {
    delete impl;
}

void pyplot::figure(const figure_options& args) const {
    auto kwargs = get_kwargs(args);
    impl->plt.attr("figure")(**kwargs);
}

void pyplot::figure(int num, const figure_options& args) const {
    auto kwargs = get_kwargs(args);
    py::object gcf = impl->pylab_helpers.attr("Gcf");
    py::object exists = gcf.attr("has_fignum")(num);
    if (exists.cast<bool>()) {
        std::cerr << "Figure " << num << " already exists\n";
        return;
    }
    impl->plt.attr("figure")(**kwargs);
}

void pyplot::subplot(int nrows, int ncols, int plot_idx) const {
    impl->plt.attr("subplot")(nrows, ncols, plot_idx);
}

void pyplot::subplot2grid(int nrows, int ncols, int rowid, int colid) const {
    py::tuple shape = py::make_tuple(nrows, ncols);
    py::tuple loc = py::make_tuple(rowid, colid);
    impl->plt.attr("subplot2grid")(shape, loc);
}

void pyplot::tight_layout() const {
    impl->plt.attr("tight_layout")();
}

void pyplot::show() const {
    impl->plt.attr("show")();
}

void pyplot::grid(bool enable) const {
    impl->plt.attr("grid")(enable);
}

void pyplot::set_aspect_equal() const {
    py::object ax = impl->plt.attr("gca")();
    ax.attr("set_aspect")("equal");
}

void pyplot::xlim(double xmin, double xmax) const {
    impl->plt.attr("xlim")(xmin, xmax);
}

void pyplot::ylim(double xmin, double xmax) const {
    impl->plt.attr("ylim")(xmin, xmax);
}

void pyplot::title(const std::string& title) const {
    impl->plt.attr("title")(title);
}

void pyplot::suptitle(const std::string& title) const {
    impl->plt.attr("suptitle")(title);
}

void pyplot::legend() const {
    impl->plt.attr("legend")();
}

void pyplot::xlabel(const std::string& text) const {
    impl->plt.attr("xlabel")(text);
}

void pyplot::ylabel(const std::string& text) const {
    impl->plt.attr("ylabel")(text);
}

void pyplot::text(double x, double y, const std::string& text) const {
    impl->plt.attr("text")(x, y, text);
}

void pyplot::marker(double x, double y, const line_options& args) const {
    plot_impl("plot", std::vector{x}, std::vector{y}, args);
}

void pyplot::arrow(double x, double y, double dx, double dy, const arrow_options& args) const {
    auto kwargs = get_kwargs(args);
    impl->plt.attr("arrow")(x, y, dx, dy, **kwargs);
}

void pyplot::plot(const array_t& x, const array_t& y, const line_options& args) const {
    plot_impl("plot", x, y, args);
}

void pyplot::plot(const array_t& x, const array_t& y, const char* format, const line_options& args) const {
    plot_impl("plot", x, y, format, args);
}

void pyplot::semilogx(const array_t& x, const array_t& y, const line_options& args) const {
    plot_impl("semilogx", x, y, args);
}

void pyplot::semilogx(const array_t& x, const array_t& y, const char* format, const line_options& args) const {
    plot_impl("semilogx", x, y, format, args);
}

void pyplot::semilogy(const array_t& x, const array_t& y, const line_options& args) const {
    plot_impl("semilogy", x, y, args);
}

void pyplot::semilogy(const array_t& x, const array_t& y, const char* format, const line_options& args) const {
    plot_impl("semilogy", x, y, format, args);
}

void pyplot::loglog(const array_t& x, const array_t& y, const line_options& args) const {
    plot_impl("loglog", x, y, args);
}

void pyplot::loglog(const array_t& x, const array_t& y, const char* format, const line_options& args) const {
    plot_impl("loglog", x, y, format, args);
}

void pyplot::fill(const array_t& x, const array_t& y, const fill_options& args) const {
    auto kwargs = get_kwargs(args);
    auto x_array = py::array_t(x.size(), x.data());
    auto y_array = py::array_t(y.size(), y.data());
    impl->plt.attr("fill")(x_array, y_array, **kwargs);
}

void pyplot::scatter(const array_t &x, const array_t &y, const scatter_options& args) const {
    auto x_array = py::array_t(x.size(), x.data());
    auto y_array = py::array_t(y.size(), y.data());
    auto kwargs = get_kwargs(args);
    impl->plt.attr("scatter")(x_array, y_array, **kwargs);
}

void pyplot::colorbar(const std::string& label) const {
    impl->plt.attr("colorbar")(py::arg("label") = label);
}

void pyplot::scatter3D(const array_t &x, const array_t &y, const array_t &z, const scatter_options& args) const {
    auto x_array = py::array_t(x.size(), x.data());
    auto y_array = py::array_t(y.size(), y.data());
    auto z_array = py::array_t(z.size(), z.data());
    auto ax = impl->axes3D();

    auto kwargs = get_kwargs(args);
    py::object mappable = ax.attr("scatter")(x_array, y_array, z_array, **kwargs);

    // Какие-то проблемные трёхмерные оси, потом найти не получается
    if (args.colorbar) {
        if (args.clabel.empty()) {
            impl->plt.attr("colorbar")(mappable, py::arg("ax") = ax);
        }
        else {
            impl->plt.attr("colorbar")(mappable, py::arg("ax") = ax, py::arg("label") = args.clabel);
        }
    }
}

void pyplot::plot_surface(const array2d_t &x, const array2d_t &y, const array2d_t &z, const surface_options& args) const {
    auto x_array = to_numpy(x);
    auto y_array = to_numpy(y);
    auto z_array = to_numpy(z);
    auto ax = impl->axes3D();
    auto kwargs = get_kwargs(args);
    ax.attr("plot_surface")(x_array, y_array, z_array, **kwargs);
}

void pyplot::plot_impl(const char* func_name, const array_t& x, const array_t& y, const line_options& args) const {
    auto x_array = py::array_t(x.size(), x.data());
    auto y_array = py::array_t(y.size(), y.data());
    auto kwargs = get_kwargs(args);
    impl->plt.attr(func_name)(x_array, y_array, **kwargs);
}

void pyplot::plot_impl(const char* func_name, const array_t& x, const array_t& y, const char* format, const line_options& args) const {
    auto x_array = py::array_t(x.size(), x.data());
    auto y_array = py::array_t(y.size(), y.data());
    auto kwargs = get_kwargs(args);
    impl->plt.attr(func_name)(x_array, y_array, format, **kwargs);
}

} // namespace zephyr::utils
#endif