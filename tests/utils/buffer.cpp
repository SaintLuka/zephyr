#include <iostream>
#include <iomanip>

#include <zephyr/utils/buffer.h>

using zephyr::utils::span;
using zephyr::utils::Buffer;
using zephyr::geom::Vector3d;

int main() {
    std::cout << std::boolalpha;

    Buffer buf = Buffer::make<double>();
    std::cout << "Buffer buf = Buffer::make<double>();\n"
        "  buf.support<double>();                // = " << buf.support<double>() << "\n" <<
        "  buf.support<double[]>();              // = " << buf.support<double[]>()  << " , can be interpreted as an array of one number\n" <<
        "  buf.support<double[1]>();             // = " << buf.support<double[1]>() << "\n" <<
        "  buf.support<double[5]>();             // = " << buf.support<double[5]>() << ", sizes mismatch\n" <<
        "  buf.support<float>();                 // = " << buf.support<float>() << "\n" <<
        "  buf.support<float[2]>();              // = " << buf.support<float[2]>() << ", checkout type is scalar\n\n";

    buf = Buffer::make<double[3]>();
    std::cout << "Buffer buf = Buffer::make<double[3]>();\n"
        "  buf.support<double>();                // = " << buf.support<double>() << "\n" <<
        "  buf.support<double[]>();              // = " << buf.support<double[]>()  << "\n" <<
        "  buf.support<double[3]>();             // = " << buf.support<double[3]>() << "\n" <<
        "  buf.support<double[5]>();             // = " << buf.support<double[5]>() << ", sizes mismatch\n" <<
        "  buf.support<Vector3d>();              // = " << buf.support<Vector3d>() << "\n" <<
        "  buf.support<std::array<double, 3>>(); // = " << buf.support<std::array<double, 3>>() << "\n\n";

    buf = Buffer::make<double[]>(3);
    std::cout << "Buffer buf = Buffer::make<double[]>(3);\n"
        "  buf.support<double>();                // = " << buf.support<double>() << "\n" <<
        "  buf.support<double[]>();              // = " << buf.support<double[]>()  << "\n" <<
        "  buf.support<double[3]>();             // = " << buf.support<double[3]>() << "\n" <<
        "  buf.support<double[5]>();             // = " << buf.support<double[5]>() << ", sizes mismatch\n" <<
        "  buf.support<Vector3d>();              // = " << buf.support<Vector3d>() << "\n" <<
        "  buf.support<std::array<double, 3>>(); // = " << buf.support<std::array<double, 3>>() << "\n";

    std::cout << "\nIf type is supported, it can be freely extracted from the buffer.\n\n";

    // При доступе к скалярному полю получаем ссылку такого типа
    buf = Buffer::make<int>(15);
    for (int i = 0; i < buf.size(); ++i) {
        buf.get_val<int>(i) = i * i;
    }
    buf.print<int>(); std::cout << "\n";

    // При доступе к векторному полю получаем span<T>
    buf = Buffer::make<int[2]>(15);
    for (int i = 0; i < buf.size(); ++i) {
        span<int> sp = buf.get_val<int[]>(i);
        sp[0] = i * i;
        sp[1] = i * i + 1;

        // Также для span есть разные плюшки вроде записи
        // из/в std::array или std::vector
    }
    buf.print<int>(); std::cout << "\n";

    // Никто не запрещает нарушать правила, поскольку внутри при доступе нет
    // никаких проверок, то можно с Vector3d работать как с double[3].

    return 0;
}
