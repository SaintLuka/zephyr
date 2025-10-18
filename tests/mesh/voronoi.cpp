// Простая демонстрация возможностей встроенного класса VDiagram,
// который реализует взвешенные диаграммы Вороного.

#include <zephyr/math/random.h>
#include <zephyr/mesh/decomp/vdiagram.h>
#include <zephyr/utils/pyplot.h>

using namespace zephyr;
using namespace zephyr::mesh;
using namespace zephyr::geom;

using decomp::VDiagram;

int main() {
    // Прямоугольная область
    Box box{Vector3d{0.0, 0.0, 0.0},
            Vector3d{2.0, 1.0, 0.0}};

    // Генератор точек в прямоугольнике
    auto rb = box.random2D(13);

    // Число ячеек
    int size = 35;

    // Генераторы ячеек диаграммы
    std::vector<Vector3d> gs(size);
    for (int i = 0; i < size; ++i) {
        gs[i] = rb.get();
    }

    // Диаграмма Вороного
    VDiagram vd(box, gs);

    // Увеличим немного параметр growth_rate, он отвечает
    // за скорость изменения весов ячеек при балансировке
    vd.growth_rate = 0.08;

    // Сбалансируем площади ячеек
    for (int i = 0; i < 10; ++i) {
        vd.balancing();
    }

    // Далее строим картинки
    utils::pyplot plt;

    plt.figure({.figsize={12.0, 7.0}});
    plt.set_aspect_equal();

    // Границы области
    plt.plot(box.outline_x(), box.outline_y(), "k");

    for (int i = 0; i < vd.size(); ++i) {
        // Границы ячеек
        plt.plot(vd.lines_x()[i], vd.lines_y()[i], "k");

        // Области сдвига генераторов (до умножения на mobility)
        plt.plot(vd.search_area_x()[i], vd.search_area_y()[i],
            {.linestyle="dotted", .linewidth=0.5, .color="black"});
    }

    // Связи между смежными ячейками диаграммы
    auto connections_x = vd.connections_x();
    auto connections_y = vd.connections_y();
    for (int i = 0; i < (int)connections_x.size(); ++i) {
        plt.plot(connections_x[i], connections_y[i],
            {.linestyle="dashed", .linewidth = 0.3, .color="blue"});
    }

    // Генераторы ячеек диаграммы
    plt.plot(vd.coords_x(), vd.coords_y(),
        {.linestyle="none", .color="tab:blue", .marker="o"});

    // Центры ячеек диаграммы
    plt.plot(vd.centers_x(), vd.centers_y(),
        {.linestyle="none", .color="tab:green", .marker="x"});

    plt.tight_layout();
    plt.show();

    return 0;
}