/// @file Тест производительности многопоточности с неравномерным доступом к памяти,
/// а также класса utils::range.
/// Здесь по логике mixed должен давать лучшее ускорение (но это не точно)

#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>

#include <zephyr/geom/vector.h>

#include <zephyr/utils/stopwatch.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/range.h>

using namespace zephyr::utils;

using Vector = zephyr::geom::Matrix<double, 15, 1>;

class Data {
public:

    Data() {
        for (int i = 0; i < ints.size(); ++i) {
            ints[i] = rand() % 20;
        }
        for (int i = 0; i < dubs.size(); ++i) {
            dubs[i] = rand() / double(RAND_MAX);
        }
    }

    std::array<int, 5> ints;
    std::array<double, 3> dubs;
};


/// @brief Представление расчетной ячейки или другого элемента
class Element {
public:
    Element() {
        res = rand() / double(RAND_MAX);
        for (int i = 0; i < a.size(); ++i) {
            a[i] = rand() / double(RAND_MAX);
            b[i] = rand() / double(RAND_MAX);
            c[i] = rand() / double(RAND_MAX);
            d[i] = rand() / double(RAND_MAX);
            e[i] = rand() / double(RAND_MAX);
            f[i] = rand() / double(RAND_MAX);
        }

        for (auto& d: temp) {
            d = rand() / double(RAND_MAX);
        }
    }

    // Некоторым образом вычисляет res
    void interact(const Data& data) {
        res = data.dubs[0] * a.dot(b) +
              data.dubs[1] * c.dot(d) +
              data.dubs[2] * e.dot(f);
    }

    std::array<double, 368> temp;
    Vector a, b, c, d, e, f;
    double res;
};

// Среднее, для проверки
double mean(std::vector<Element>& elements) {
    double res = 0.0;
    for (auto &elem: elements) {
        res += elem.res;
    }
    return res / elements.size();
}

// Среднее, для проверки
double mean(std::vector<std::pair<Element, Data>>& elements) {
    double res = 0.0;
    for (auto &elem: elements) {
        res += elem.first.res;
    }
    return res / elements.size();
}

void foo_range(size_t idx, std::vector<Element>& elems, std::vector<Data>& datas) {
    elems[idx].interact(datas[idx]);
}

void foo_mixed(std::pair<Element, Data>& p) {
    p.first.interact(p.second);
}


int main() {
    int n_elements = 100000;
    int n_steps = 500;

    std::vector<Element> elements_1(n_elements);
    std::vector<Data> data_1(n_elements);

    std::vector<std::pair<Element, Data>> mixed_1(n_elements);
    for (int i = 0; i < n_elements; ++i) {
        mixed_1[i].first = elements_1[i];
        mixed_1[i].second = data_1[i];
    }

    std::vector<Data> data_2(data_1);
    std::vector<Element> elements_2(elements_1);
    std::vector<std::pair<Element, Data>> mixed_2(mixed_1);

    std::cout << "    elems_1.mean: " << mean(elements_1) << "\n";
    std::cout << "    mixed_1.mean: " << mean(mixed_1) << "\n";
    std::cout << "    elems_2.mean: " << mean(elements_2) << "\n";
    std::cout << "    mixed_2.mean: " << mean(mixed_2) << "\n\n";

    // Последовательный запуск
    
    Stopwatch sw_range_single(true);
    for (int step: range(n_steps)) {
        threads::for_each(range(n_elements).begin(), range(n_elements).end(), foo_range,
                          std::ref(elements_1), std::ref(data_1));
    }
    sw_range_single.stop();
    
    Stopwatch sw_mixed_single(true);
    for (int step: range(n_steps)) {
        threads::for_each(mixed_1.begin(), mixed_1.end(), foo_mixed);
    }
    sw_mixed_single.stop();

    std::cout << "    elems_1.mean: " << mean(elements_1) << "\n";
    std::cout << "    mixed_1.mean: " << mean(mixed_1) << "\n";

    // Многопоточный запуск

    threads::on();

    Stopwatch sw_range_thread(true);
    for (int step: range(n_steps)) {
        threads::for_each(range(n_elements).begin(), range(n_elements).end(), foo_range,
                          std::ref(elements_2), std::ref(data_2));
    }
    sw_range_thread.stop();

    Stopwatch sw_mixed_thread(true);
    for (int step: range(n_steps)) {
        threads::for_each(mixed_2.begin(), mixed_2.end(), foo_mixed);
    }
    sw_mixed_thread.stop();

    std::cout << "    elems_2.mean: " << mean(elements_2) << "\n";
    std::cout << "    mixed_2.mean: " << mean(mixed_2) << "\n\n";


    // Вывод статистики

    std::cout << std::setprecision(3);

    std::cout << "          |  Serial, ms  | Parallel, ms | Speed Up |\n";
    std::cout << "----------+--------------+--------------+----------|\n";

    std::cout << "  range   | "
              << std::setw(12) << sw_range_single.milliseconds() << " | "
              << std::setw(12) << sw_range_thread.milliseconds() << " | "
              << std::setw(8) << double(sw_range_single.milliseconds()) / double(sw_range_thread.milliseconds()) << " |\n";

    std::cout << "  mixed   | "
              << std::setw(12) << sw_mixed_single.milliseconds() << " | "
              << std::setw(12) << sw_mixed_thread.milliseconds() << " | "
              << std::setw(8) << double(sw_mixed_single.milliseconds()) / double(sw_mixed_thread.milliseconds()) << " |\n";

    return 0;
}