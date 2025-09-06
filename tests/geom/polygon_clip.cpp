// Тестирование класса Polygon. Проверка функций интегрирования по полигону,
// отсечения от полигона линией и кругом.

#include <iostream>

#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/triangle.h>

#include <zephyr/utils/matplotlib.h>

using namespace zephyr::geom;
namespace plt = zephyr::utils::matplotlib;


// Тестовая прямая
struct TLine {
    Vector3d p;  //< Точка прямой
    Vector3d n;  //< Внешняя нормаль

    TLine(const Vector3d& _p, const Vector3d& _n)
            : p(_p), n(_n) { }

    // характеристическая функция
    bool operator()(const Vector3d& v) const {
        return (v - p).dot(n) < 0.0;
    };
};

// Тестовая окружность
struct TDisk {
    Vector3d c;  //< Центр окружности
    double r;    //< Радиус окружности

    TDisk(const Vector3d& _c, double _r)
        : c(_c), r(_r) { }

    // характеристическая функция
    bool operator()(const Vector3d& p) const {
        return (p - c).squaredNorm() < r * r;
    };
};

int main() {
    // Тестовый многоугольник
    Polygon poly = {
            Vector3d{-5.0, -3.0, 0.0},
            Vector3d{+6.0, -4.0, 0.0},
            Vector3d{+5.0, +4.0, 0.0},
            Vector3d{-2.0, +4.0, 0.0},
            Vector3d{-7.0, +1.0, 0.0}
    };
    poly.sort();

    std::cout << "Polygon area:  " << poly.area() << "\n\n";

    // Сечение многоугольника прямой
    if (true) {
        // Тестовая прямая
        TLine line(poly[0] + 0.28 * (poly[2] - poly[0]),
                   (poly[2] - poly[0]).normalized());

        Polygon clip = poly.clip(line.p, line.n);

        double S1 = poly.clip_area(line);
        double S2 = clip.area();
        double S3 = poly.clip_area(line.p, line.n);

        // Находим объемную долю сечения
        double vf = S3 / poly.area();

        // Восстанавливаем сечение по объемной доле
        Vector3d p = poly.find_section(line.n, vf);
        double S4 = poly.clip_area(p, line.n);

        std::cout << "  Line clip S1: " << S1 << ";\terr: " << std::abs(S1 - S3) << "\t(approx)\n";
        std::cout << "  Line clip S2: " << S2 << ";\terr: " << std::abs(S2 - S3) << "\n";
        std::cout << "  Line clip S3: " << S3 << ";\terr: " << std::abs(S3 - S3) << "\n";
        std::cout << "  Line clip S4: " << S4 << ";\terr: " << std::abs(S4 - S3) << "\n\n";

        std::cout << "  Volume fraction: " << vf << "\n";
        std::cout << "  Section point: " << p.transpose() << ";\terr: " << (p - line.p).dot(line.n) << "\n\n";

        // Строим многоугольник и сечение
        plt::figure_size(14.0, 5.0);
        plt::subplot(1, 2, 1);
        plt::title("Сечение многоугольника прямой линией");
        plt::set_aspect_equal();
        plt::plot(poly.xs(), poly.ys());
        plt::fill(clip.xs(), clip.ys(), {{"color", "#0000ff3f"}});
    }

    // Сечение многоугольника окружностью
    if (true) {
        // Тестовая окружность
        TDisk disk({1.5, -0.3, 0.0}, 5.0);

        // Построение приближенного сечения
        int M = 10000;
        std::vector<Vector3d> circle;
        circle.reserve(M + poly.size());
        for (int i = 0; i < M; ++i) {
            double phi = 2 * M_PI * i / (M - 1.0);
            Vector3d n = { std::cos(phi), std::sin(phi), 0.0 };
            Vector3d p = disk.c + disk.r * n;

            if (poly.inside(p)) {
                circle.push_back(p);
            }
        }
        for (int i = 0; i < poly.size(); ++i) {
            if (disk(poly[i])) {
                circle.push_back(poly[i]);
            }
        }
        Polygon dclip(circle);
        dclip.sort();

        double S1 = poly.clip_area(disk);
        double S2 = dclip.area();
        double S3 = poly.disk_clip_area(disk.c, disk.r);

        // Находим объемную долю сечения
        double vf = S3 / poly.area();

        // Найдем нормаль и соответствующее сечение прямой
        Vector3d n = poly.disk_clip_normal(disk.c, disk.r);
        Vector3d p = poly.find_section(n, vf);
        Polygon clip = poly.clip(p, n);
        double S4 = poly.clip_area(p, n);
        double S5 = clip.area();

        std::cout << "  Disk clip S1: " << S1 << ";\terr: " << std::abs(S1 - S3) << "\t(approx)\n";
        std::cout << "  Disk clip S2: " << S2 << ";\terr: " << std::abs(S2 - S3) << "\t(approx)\n";
        std::cout << "  Disk clip S3: " << S3 << ";\terr: " << std::abs(S3 - S3) << "\n";
        std::cout << "  Disk clip S4: " << S4 << ";\terr: " << std::abs(S4 - S3) << "\n";
        std::cout << "  Disk clip S5: " << S5 << ";\terr: " << std::abs(S5 - S3) << "\n\n";

        std::cout << "  Volume fraction: " << vf << "\n";
        std::cout << "  Section point: " << p.transpose() << "\n\n";

        // Строим многоугольник и сечение
        plt::subplot(1, 2, 2);
        plt::title("Сечение многоугольника окружностью");
        plt::set_aspect_equal();
        plt::plot(poly.xs(), poly.ys());
        //plt::plot({disk.c.x()}, {disk.c.y()}, "kx");
        plt::fill(dclip.xs(), dclip.ys(), {{"color", "#00ff003f"}});
        plt::fill(clip.xs(), clip.ys(), {{"color", "#0000ff3f"}});
    }

    plt::tight_layout();
    plt::show();

    return 0;
}