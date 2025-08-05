// @brief Тест производительности многопоточности при работе с векторными
// контейнерами трех типов. Проверяются следующие типы:
// std::array -- Классический массив с фиксированным размером.
// std::vector -- Классический динамический массив из стандартной библиотеки.
// boost::container::static_vector -- Вектор с ограниченным максимальным
// размером, известным на момент компиляции. Использует буфер памяти
// константного размера. Имитирует функционал std::vector.
// boost::container::small_vector -- Вектор содержит буфер константного
// размера, при превышении размеров данные переносятся в динамическую память.
//
// Задача возникла при расчетах многоматериальных течений, в данных задачах
// создается множество векторов одинакового размера, с которыми производятся
// арифметические операции. При этом точный размер векторов не известен на
// этапе компиляции.
//
// Надо тестировать на разных платформах, но предварительные выводы такие:
// small_vector -- совсем обсосный
// static_vector лучше std::vector даже в однопоточном режиме, особенно
// при размерах вектора около 5 double (выигрывает примерно в 10 раз),
// но при размерах в 50 double различия не заметны.
//
// На домашнем (8 ядер) вообще никакой разницы, даже std::vector побеждает
// при размерах в 15 double. Ничего не понятно, std::vector очень мощный,
// может зря я от него отказывался.
//
// Добавил std::array, теперь он побеждает, фигли он быстрее static_vector?
// По логике должно быть совершенно одно и то же. Надо ещё vector из Eigen
// добавить, да чтож они все по-разному работают??

#include <iostream>
#include <iomanip>
#include <cmath>

#include <array>
#include <functional>
#include <vector>
#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>

#include <zephyr/utils/stopwatch.h>
#include <zephyr/utils/threads.h>

using namespace zephyr::utils;

// Актуальные размеры: 15 - 25 для PT-модели
// ~ 15 * n_materials для полностью неравновесных
static constexpr int SIZE = 15;

using std_array = std::array<double, SIZE>;
using std_vector = std::vector<double>;
using small_vector = boost::container::small_vector<double, SIZE>;
using static_vector = boost::container::static_vector<double, SIZE>;

template <typename container>
void resize(container& arr) {
    arr.resize(SIZE);
}

template <>
void resize<std_array>(std_array& arr) {

}

// Простой класс для хранения вектора
// container Тип хранения вектора
template <typename container>
class Vector {
public:
    Vector() {
        resize(vec);
    }

    Vector(const Vector& v)
        : vec(v.vec) { }

    static Vector Zero() {
        return Vector();
    }

    static Vector Rand() {
        Vector vec;
        for (int i = 0; i < SIZE; ++i) {
            vec[i] = 2.0 * rand() / double(RAND_MAX) - 1.0;
        }
        return vec;
    }

    int size() const {
        return vec.size();
    }

    double& operator[](int idx) {
        return vec[idx];
    }

    const double& operator[](int idx) const {
        return vec[idx];
    }

    // Бессмысленная операция, проекция
    Vector dot(const Vector& v) const {
        double coeff = 0.0;
        for (int i = 0; i < SIZE; ++i) {
            coeff += vec[i] * v[i];
        }
        Vector res;
        for (int i = 0; i < SIZE; ++i) {
            res[i] = vec[i] * coeff;
        }
        return res;
    }

    // Сложение
    Vector operator+(const Vector& v) {
        Vector res = Vector::Zero();
        for (int i = 0; i < size(); ++i) {
            res[i] = vec[i] + v.vec[i];
        }
        return res;
    }

    void operator+=(const Vector& v) {
        for (int i = 0; i < size(); ++i) {
            vec[i] += v.vec[i];
        }
    }

    // Умножение на число
    Vector operator*(double q) {
        Vector res = Vector::Zero();
        for (int i = 0; i < size(); ++i) {
            res[i] = vec[i] * q;
        }
        return res;
    }
    
    double mean() const {
        double res = 0.0;
        for (int i = 0; i < size(); ++i) {
            res += vec[i];
        }
        return res / SIZE;
    }

protected:
    container vec;
};

template <typename container>
std::ostream& operator<<(std::ostream& os, const Vector<container>& vec) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "{ ";
    for (int i = 0; i < vec.size() - 1; ++i) {
        os << std::setw(6) << vec[i] << ", ";
    }
    os << vec[vec.size() - 1] << " }";
    return os;
}

using StdArray = Vector<std_array>;
using StdVector = Vector<std_vector>;
using SmallVector = Vector<small_vector>;
using StaticVector = Vector<static_vector>;


// Представление расчетной ячейки или другого элемента
// Vec Тип вектора
template <typename Vec>
class Element {
public:
    Element() {
        a = Vec::Rand();
        b = Vec::Rand();
        c = Vec::Rand();
        d = Vec::Rand();
        e = Vec::Rand();
        f = Vec::Rand();
    }
    
    double mean() const {
        double res = 0.0;
        res += a.mean();
        res += b.mean();
        res += c.mean();
        res += d.mean();
        res += e.mean();
        res += f.mean();
        return res / 6.0;
    }
    
    std::array<double, 368> data;
    Vec a, b, c, d, e, f;
};

using StdElement = Element<Vector<std_vector>>;
using SmallElement = Element<Vector<small_vector>>;
using StaticElement = Element<Vector<static_vector>>;
using ArrayElement = Element<Vector<std_array>>;

// Среднее, для проверки
template <typename E>
double mean(std::vector<E>& elements) {
    double res = 0.0;
    for (auto &elem: elements) {
        res += elem.mean();
    }
    return res / elements.size();
}

// Функция 1. Обработать один элемент независимо от остальных
template <typename Vec>
void foo_1(Element<Vec>& elem) {
    Vec a = elem.a + elem.b * 1.1;
    Vec b = elem.c + elem.d * 2.3;
    Vec c = elem.e + elem.f * 4.5;

    elem.data[0] = a.dot(b * 0.1).dot(c * 0.1).mean();
    elem.data[1] = b.dot(c * 0.2).dot(a * 0.2).mean();
    elem.data[2] = c.dot(a * 0.3).dot(b * 0.3).mean();
}

// Функция 2. Обработать один элемент на основе пары других элементов
// из массива. Моделирует взаимодействие элементов.
template <typename Vec>
void foo_2(Element<Vec>& elem, std::vector<Element<Vec>>& elems) {
    int ic = &elem - &elems[0];
    int shift = 1; //elems.size() / 3;
    int ir = (ic + shift + elems.size()) % elems.size();
    int il = (ic - shift + elems.size()) % elems.size();

    Element<Vec> &elem_L = elems[il];
    Element<Vec> &elem_R = elems[ir];

    elem.data[0] += (elem_L.a + elem_R.a).dot(elem.a * 0.1).mean();
    elem.data[1] += (elem_L.b + elem_R.b).dot(elem.b * 0.2).mean();
    elem.data[2] += (elem_L.c + elem_R.c).dot(elem.c * 0.3).mean();
    elem.data[3] += (elem_L.d + elem_R.d).dot(elem.d * 0.4).mean();
    elem.data[4] += (elem_L.e + elem_R.e).dot(elem.e * 0.5).mean();
    elem.data[5] += (elem_L.f + elem_R.f).dot(elem.f * 0.6).mean();
}


int main() {
    std::cout << "Generate three arrays with different container types.\n";
    std::cout << "  Elements are filled random, checkout mean.\n";
    
    int n_elements = 200000;
    int n_steps = 50;
    
    srand(-13);
    std::vector<StdElement> elements_1(n_elements);

    srand(-13);
    std::vector<SmallElement> elements_2(n_elements);

    srand(-13);
    std::vector<StaticElement> elements_3(n_elements);

    srand(-13);
    std::vector<ArrayElement> elements_4(n_elements);

    std::cout << "    elements_1.mean: " << mean(elements_1) << "\n";
    std::cout << "    elements_2.mean: " << mean(elements_2) << "\n";
    std::cout << "    elements_3.mean: " << mean(elements_3) << "\n";
    std::cout << "    elements_4.mean: " << mean(elements_4) << "\n\n";

    // Последовательный запуск
    
    Stopwatch sw_std_single_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_1.begin(), elements_1.end(), foo_1<StdVector>);
    }
    sw_std_single_1.stop();
    
    Stopwatch sw_small_single_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_2.begin(), elements_2.end(), foo_1<SmallVector>);
    }
    sw_small_single_1.stop();

    Stopwatch sw_static_single_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_3.begin(), elements_3.end(), foo_1<StaticVector>);
    }
    sw_static_single_1.stop();

    Stopwatch sw_array_single_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_4.begin(), elements_4.end(), foo_1<StdArray>);
    }
    sw_array_single_1.stop();
    
    
    Stopwatch sw_std_single_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_1.begin(), elements_1.end(), foo_2<StdVector>, std::ref(elements_1));
    }
    sw_std_single_2.stop();

    Stopwatch sw_small_single_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_2.begin(), elements_2.end(), foo_2<SmallVector>, std::ref(elements_2));
    }
    sw_small_single_2.stop();

    Stopwatch sw_static_single_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_3.begin(), elements_3.end(), foo_2<StaticVector>, std::ref(elements_3));
    }
    sw_static_single_2.stop();

    Stopwatch sw_array_single_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_4.begin(), elements_4.end(), foo_2<StdArray>, std::ref(elements_4));
    }
    sw_array_single_2.stop();

    std::cout << "  After first operation:\n";
    std::cout << "    elements_1.mean: " << mean(elements_1) << "\n";
    std::cout << "    elements_2.mean: " << mean(elements_2) << "\n";
    std::cout << "    elements_3.mean: " << mean(elements_3) << "\n";
    std::cout << "    elements_4.mean: " << mean(elements_4) << "\n\n";


    // Многопоточный запуск

    threads::on();

    Stopwatch sw_std_thread_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_1.begin(), elements_1.end(), foo_1<StdVector>);
    }
    sw_std_thread_1.stop();

    Stopwatch sw_small_thread_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_2.begin(), elements_2.end(), foo_1<SmallVector>);
    }
    sw_small_thread_1.stop();

    Stopwatch sw_static_thread_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_3.begin(), elements_3.end(), foo_1<StaticVector>);
    }
    sw_static_thread_1.stop();

    Stopwatch sw_array_thread_1(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_4.begin(), elements_4.end(), foo_1<StdArray>);
    }
    sw_array_thread_1.stop();


    Stopwatch sw_std_thread_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_1.begin(), elements_1.end(), foo_2<StdVector>, std::ref(elements_1));
    }
    sw_std_thread_2.stop();

    Stopwatch sw_small_thread_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_2.begin(), elements_2.end(), foo_2<SmallVector>, std::ref(elements_2));
    }
    sw_small_thread_2.stop();

    Stopwatch sw_static_thread_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_3.begin(), elements_3.end(), foo_2<StaticVector>, std::ref(elements_3));
    }
    sw_static_thread_2.stop();

    Stopwatch sw_array_thread_2(true);
    for(int step = 0; step < n_steps; ++step) {
        threads::for_each(elements_4.begin(), elements_4.end(), foo_2<StdArray>, std::ref(elements_4));
    }
    sw_array_thread_2.stop();

    std::cout << "  After second operation:\n";
    std::cout << "    elements_1.mean: " << mean(elements_1) << "\n";
    std::cout << "    elements_2.mean: " << mean(elements_2) << "\n";
    std::cout << "    elements_3.mean: " << mean(elements_3) << "\n";
    std::cout << "    elements_4.mean: " << mean(elements_4) << "\n\n";


    // Вывод статистики

    std::cout << std::setprecision(3);

    std::cout << "                |  Serial, ms  | Parallel, ms | Speed Up |\n";
    std::cout << "----------------+--------------+--------------+----------|\n";

    std::cout << "foo_1:\n";
    std::cout << "  std::vector   | "
              << std::setw(12) << sw_std_single_1.milliseconds() << " | "
              << std::setw(12) << sw_std_thread_1.milliseconds() << " | "
              << std::setw(8) << double(sw_std_single_1.milliseconds()) / double(sw_std_thread_1.milliseconds()) << " |\n";

    std::cout << "  small_vector  | "
              << std::setw(12) << sw_small_single_1.milliseconds() << " | "
              << std::setw(12) << sw_small_thread_1.milliseconds() << " | "
              << std::setw(8) << double(sw_small_single_1.milliseconds()) / double(sw_small_thread_1.milliseconds()) << " |\n";

    std::cout << "  static_vector | "
              << std::setw(12) << sw_static_single_1.milliseconds() << " | "
              << std::setw(12) << sw_static_thread_1.milliseconds() << " | "
              << std::setw(8) << double(sw_static_single_1.milliseconds()) / double(sw_static_thread_1.milliseconds()) << " |\n";

    std::cout << "  std::array    | "
              << std::setw(12) << sw_array_single_1.milliseconds() << " | "
              << std::setw(12) << sw_array_thread_1.milliseconds() << " | "
              << std::setw(8) << double(sw_array_single_1.milliseconds()) / double(sw_array_thread_1.milliseconds()) << " |\n";

    std::cout << "\nfoo_2:\n";
    std::cout << "  std::vector   | "
              << std::setw(12) << sw_std_single_2.milliseconds() << " | "
              << std::setw(12) << sw_std_thread_2.milliseconds() << " | "
              << std::setw(8) << double(sw_std_single_2.milliseconds()) / double(sw_std_thread_2.milliseconds()) << " |\n";

    std::cout << "  small_vector  | "
              << std::setw(12) << sw_small_single_2.milliseconds() << " | "
              << std::setw(12) << sw_small_thread_2.milliseconds() << " | "
              << std::setw(8) << double(sw_small_single_2.milliseconds()) / double(sw_small_thread_2.milliseconds()) << " |\n";

    std::cout << "  static_vector | "
              << std::setw(12) << sw_static_single_2.milliseconds() << " | "
              << std::setw(12) << sw_static_thread_2.milliseconds() << " | "
              << std::setw(8) << double(sw_static_single_2.milliseconds()) / double(sw_static_thread_2.milliseconds()) << " |\n";

    std::cout << "  std::array    | "
              << std::setw(12) << sw_array_single_2.milliseconds() << " | "
              << std::setw(12) << sw_array_thread_2.milliseconds() << " | "
              << std::setw(8) << double(sw_array_single_2.milliseconds()) / double(sw_array_thread_2.milliseconds()) << " |\n";

    return 0;
}