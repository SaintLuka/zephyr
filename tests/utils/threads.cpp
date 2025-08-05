// @brief Бенчмарк для тестирования производительности в многопоточном режиме.
// Производительность тестируется при обработке массива, который содержит
// матрицы относительно большого размера, в качестве трудоемкой операции
// используется поиск собственных значений.

#include <iostream>
#include <iomanip>
#include <cmath>

#include <Eigen/Dense>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::utils;

// @return Случайное число от [-1, 1]
inline double urand() {
    return 2.0 * double(rand()) / RAND_MAX - 1.0;
}

// @return Функция знака
inline double sign(double x) {
    return x < 0.0 ? -1.0 : 1.0;
}

const int size = 20; //< Размерность используемых матриц

using Vector = Eigen::Matrix<double, 1, size>;
using Matrix = Eigen::Matrix<double, size, size>;

// @brief Тип для хранения
struct Element {
    Matrix mat;

    // @brief Создать элемент по матрице
    explicit Element(const Matrix& _mat) : mat(_mat) { }

    // @brief Создать случайную матрицу size x size
    // с вещественными собственными значениями от -1 до 1.
    Element() {
        Matrix A;
        Vector L;

        for (int i = 0; i < size; ++i) {
            L(i) = urand();
            for (int j = 0; j < size; ++j) {
                A(i, j) = urand();
            }
        }

        mat = A.inverse() * L.asDiagonal() * A;
    }

    // @brief Трудоемкая операция поворота,
    // собственные значения матрицы не изменяются
    void rotate() {
        Vector L = mat.eigenvalues().real();
        Matrix R = Vector::Ones().asDiagonal();

        // Генерируем матрицу вращения
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                double phi = 0.5 * M_PI * (L(i) + L(j));

                // Матрица вращения в плоскости осей (i, j)
                Matrix A = Vector::Ones().asDiagonal();

                A(i, i) = std::cos(phi);
                A(j, j) = std::cos(phi);
                A(i, j) = std::sin(phi);
                A(j, i) = -std::sin(phi);

                R *= A;
            }
        }

        // Вращение матрицы
        mat = R * mat * R.transpose();
    }

    // @return Среднее собственное значение
    double avg_eigen() const {
        return mat.eigenvalues().mean().real();
    }

    // @brief Оператор сравнения
    bool operator<(const Element& r) const {
        return avg_eigen() < r.avg_eigen();
    }

    // @brief Оператор сравнения
    bool operator>(const Element& r) const {
        return avg_eigen() > r.avg_eigen();
    }

    // @brief Корень степени K от каждого элемента матрицы
    static Matrix reduce(Matrix& m, int K) {
        Matrix res;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                res(i, j) = sign(m(i, j)) * std::pow(std::abs(m(i, j)), 1.0 / K);
            }
        }
        return res;
    }

    // @brief Операция суммирования
    Element& operator+=(const Element& r) {
        mat += r.mat + 0.0 * Matrix(r.mat.eigenvalues().real().asDiagonal());
        return *this;
    }

    // @brief Операция свертки (покомпонентный минимум)
    Element& operator&=(const Element& r) {
        mat = mat.cwiseMin(r.mat) + 0.0 * Matrix(r.mat.eigenvalues().real().asDiagonal());
        return *this;
    }

    // @brief Инициализация при поиске минимума
    static Element init_min() {
        return Element(+100.0 * Vector::Ones().asDiagonal());
    }

    // @brief Инициализация при поиске максимума
    static Element init_max() {
        return Element(-100.0 * Vector::Ones().asDiagonal());
    }

    // @brief Инициализация при суммировании
    static Element init_sum() {
        return Element(Matrix::Zero());
    }

    // @brief Инициализация при свертке
    static Element init_reduce() {
        return Element(1.0e300 * Matrix::Ones());
    }

    // @brief Вытащить блок 3 х 3
    Eigen::Matrix3d block() {
        return Eigen::Matrix3d(mat.block<3, 3>(0, 0));
    }
};

struct Info {
    Stopwatch sw_full;
    Stopwatch sw_rotate;
    Stopwatch sw_min;
    Stopwatch sw_max;
    Stopwatch sw_sum;
    Stopwatch sw_reduce;

    Element elem_min;
    Element elem_max;
    Element elem_sum;
    Element elem_reduce;

    Info() = default;
};

Info foo1(std::vector<Element>& elements) {
    Info info;

    info.sw_full.start();

    // Поворачиваем все элементы
    info.sw_rotate.start();
    for (auto& elem: elements) {
        elem.rotate();
    }
    info.sw_rotate.stop();

    // Находим минимальный
    info.sw_min.start();
    info.elem_min = Element::init_min();
    for (auto& elem: elements) {
        if (elem < info.elem_min) {
            info.elem_min = elem;
        }
    }
    info.sw_min.stop();

    // Находим максимальный
    info.sw_max.start();
    info.elem_max = Element::init_max();
    for (auto& elem: elements) {
        if (elem > info.elem_max) {
            info.elem_max = elem;
        }
    }
    info.sw_max.stop();

    // Находим сумму
    info.sw_sum.start();
    info.elem_sum = Element::init_sum();
    for (auto& elem: elements) {
        info.elem_sum += elem;
    }
    info.elem_sum.mat /= elements.size();
    info.sw_sum.stop();

    // Находим свертку
    info.sw_reduce.start();
    info.elem_reduce = Element::init_reduce();
    for (auto& elem: elements) {
        info.elem_reduce &= elem;
    }
    info.sw_reduce.stop();

    info.sw_full.stop();

    return info;
}

Info foo2(std::vector<Element>& elements) {
    using Iterator = std::vector<Element>::iterator;

    // Возвращает сам элемент по интератору
    auto get_elem = [](Element& elem) -> Element {
        return elem;
    };

    Info info;

    info.sw_full.start();

    // Поворачиваем все элементы
    info.sw_rotate.start();
    threads::for_each(
            elements.begin(), elements.end(),
            [](Element& elem) {
                elem.rotate();
            });
    info.sw_rotate.stop();

    // Находим минимальный
    info.sw_min.start();
    info.elem_min = threads::min(
            elements.begin(), elements.end(),
            Element::init_min(), get_elem);
    info.sw_min.stop();

    // Находим максимальный
    info.sw_max.start();
    info.elem_max = threads::max(
            elements.begin(), elements.end(),
            Element::init_max(), get_elem);
    info.sw_max.stop();

    // Находим сумму
    info.sw_sum.start();
    info.elem_sum = threads::sum(
            elements.begin(), elements.end(),
            Element::init_sum(), get_elem);
    info.elem_sum.mat /= elements.size();
    info.sw_sum.stop();

    // Находим свертку
    info.sw_reduce.start();
    info.elem_reduce = threads::reduce(
            elements.begin(), elements.end(),
            Element::init_reduce(), get_elem);
    info.sw_reduce.stop();

    info.sw_full.stop();

    return info;
}

int main() {
    int n_elements = 10000;

    std::vector<Element> arr1(n_elements);
    std::vector<Element> arr2(arr1);
    std::vector<Element> arr3(arr1);

    // Вычисляется последовательно
    Info info1 = foo1(arr1);

    // Вычисляется последовательно, но через интерфейс threads
    Info info2 = foo2(arr2);

    // Треды включены, вычисляется в много потоков
    threads::on();
    Info info3 = foo2(arr3);


    // Вывод результатов

    std::cout << "arr[0]:\n";
    std::cout << arr1[0].block() << "\n\n";

    std::cout << "arr.min:\n";
    std::cout << " avg_eigen: " << info1.elem_min.avg_eigen() << "\n";
    std::cout << info1.elem_min.block() << "\n\n";

    std::cout << "arr.max:\n";
    std::cout << " avg_eigen: " << info1.elem_max.avg_eigen() << "\n";
    std::cout << info1.elem_max.block() << "\n\n";

    std::cout << "arr.sum:\n";
    std::cout << info1.elem_sum.block() << "\n\n";

    std::cout << "arr.reduce:\n";
    std::cout << info1.elem_reduce.block() << "\n\n";


    // Проверяем, что результаты совпадают

    std::cout << "n_threads: " << threads::count() << "\n\n";

    // Матрицы после поворотов совпадают
    {
        double rot_err1(0.0);
        double rot_err2(0.0);
        for (int i = 0; i < n_elements; ++i) {
            rot_err1 = std::max(rot_err1, (arr2[i].mat - arr1[i].mat).cwiseAbs().maxCoeff());
            rot_err2 = std::max(rot_err2, (arr3[i].mat - arr1[i].mat).cwiseAbs().maxCoeff());
        }
        std::cout << "Rotation  operation error:  " << rot_err1 << " " << rot_err2 << "\n";
    }

    // Минимум совпадает
    {
        double min_err1 = (info2.elem_min.mat - info1.elem_min.mat).cwiseAbs().maxCoeff();
        double min_err2 = (info3.elem_min.mat - info1.elem_min.mat).cwiseAbs().maxCoeff();
        std::cout << "Minimum   operation error:  " << min_err1 << " " << min_err2 << "\n";
    }

    // Максимум совпадает
    {
        double max_err1 = (info2.elem_max.mat - info1.elem_max.mat).cwiseAbs().maxCoeff();
        double max_err2 = (info3.elem_max.mat - info1.elem_max.mat).cwiseAbs().maxCoeff();
        std::cout << "Maximum   operation error:  " << max_err1 << " " << max_err2 << "\n";
    }

    // Сумма совпадает
    {
        double sum_err1 = (info2.elem_sum.mat - info1.elem_sum.mat).cwiseAbs().maxCoeff();
        double sum_err2 = (info3.elem_sum.mat - info1.elem_sum.mat).cwiseAbs().maxCoeff();
        std::cout << "Summation operation error:  " << sum_err1 << " " << sum_err2 << "\n";
    }

    // Свертка совпадает
    {
        double red_err1 = (info2.elem_reduce.mat - info1.elem_reduce.mat).cwiseAbs().maxCoeff();
        double red_err2 = (info3.elem_reduce.mat - info1.elem_reduce.mat).cwiseAbs().maxCoeff();
        std::cout << "Reduce    operation error:  " << red_err1 << " " << red_err2 << "\n\n";
    }


    // Таблица производительности.
    // Строка elapsed - Полное время выполнения
    // Столбец Serial(1) - Последовательное выполенени
    // Столбец Serial(2) - Последовательное выполнение через интерфейс тредов
    // Столбец Parallel - Параллельное выполнение на тредах
    // Столбцы Speed Up - Ускорение выполнения по сравнению со столбцом Serial(2)
    // Speed Up для Serial(2) должен быть около единицы, но не сильно меньше.

    std::cout << std::setprecision(3);

    std::cout << "            | Serial(1), ms | Serial(2), ms | Speed Up | Parallel, ms | Speed Up |\n";
    std::cout << "------------+---------------+---------------+----------+--------------+----------|\n";
    
    std::cout << "Elapsed     | "
              << std::setw(13) << info1.sw_full.milliseconds() << " | "
              << std::setw(13) << info2.sw_full.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_full.milliseconds()) / info2.sw_full.milliseconds() << " | "
              << std::setw(12) << info3.sw_full.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_full.milliseconds()) / info3.sw_full.milliseconds() << " |\n";
    
    std::cout << "  foreach   | "
              << std::setw(13) << info1.sw_rotate.milliseconds() << " | "
              << std::setw(13) << info2.sw_rotate.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_rotate.milliseconds()) / info2.sw_rotate.milliseconds() << " | "
              << std::setw(12) << info3.sw_rotate.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_rotate.milliseconds()) / info3.sw_rotate.milliseconds() << " |\n";
    
    std::cout << "  minimum   | "
              << std::setw(13) << info1.sw_min.milliseconds() << " | "
              << std::setw(13) << info2.sw_min.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_min.milliseconds()) / info2.sw_min.milliseconds() << " | "
              << std::setw(12) << info3.sw_min.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_min.milliseconds()) / info3.sw_min.milliseconds() << " |\n";
    
    std::cout << "  maximum   | "
              << std::setw(13) << info1.sw_max.milliseconds() << " | "
              << std::setw(13) << info2.sw_max.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_max.milliseconds()) / info2.sw_max.milliseconds() << " | "
              << std::setw(12) << info3.sw_max.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_max.milliseconds()) / info3.sw_max.milliseconds() << " |\n";
    
    std::cout << "  summation | "
              << std::setw(13) << info1.sw_sum.milliseconds() << " | "
              << std::setw(13) << info2.sw_sum.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_sum.milliseconds()) / info2.sw_sum.milliseconds() << " | "
              << std::setw(12) << info3.sw_sum.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_sum.milliseconds()) / info3.sw_sum.milliseconds() << " |\n";
    
    std::cout << "  reduce    | "
              << std::setw(13) << info1.sw_reduce.milliseconds() << " | "
              << std::setw(13) << info2.sw_reduce.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_reduce.milliseconds()) / info2.sw_reduce.milliseconds() << " | "
              << std::setw(12) << info3.sw_reduce.milliseconds() << " | "
              << std::setw(8) << double(info1.sw_reduce.milliseconds()) / info3.sw_reduce.milliseconds() << " |\n";

    return 0;
}