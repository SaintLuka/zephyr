#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/euler/soa_mesh.h>

using zephyr::geom::Vector3d;
using namespace zephyr::mesh;


// Имитация вектора состояния
struct PState {
    Storable<double>   rho;  // Плотность
    Storable<Vector3d> v;    // Скорость
    Storable<double>   P;    // Давление
    Storable<double>   e;    // Энергия

    // Изменяет хранилище s, добавляет в него новые поля
    void add_types(SoaStorage &s, const std::string &prefix) {
        rho = s.add<double>  (prefix + ".rho");
        v   = s.add<Vector3d>(prefix + ".v");
        P   = s.add<double>  (prefix + ".P");
        e   = s.add<double>  (prefix + ".e");
    }
};

// Имитация расширенного вектора состояния
struct State {
    PState curr;  // Текущее состояние
    PState half;  // Состояние на полушаге
    PState next;  // Состояние на следующем шаге

    /// @brief Градиент вектора состояния
    PState d_dx, d_dy, d_dz;

    // Изменяет хранилище s, добавляет в него новые поля
    void add_types(SoaStorage &s) {
        curr.add_types(s, "curr");
        half.add_types(s, "half");
        next.add_types(s, "next");
        d_dx.add_types(s, "d_dx");
        d_dy.add_types(s, "d_dy");
        d_dz.add_types(s, "d_dz");
    }
};

inline std::string yes_no(bool val) { return val ? "yes" : "no"; }


int main() {
    SoaStorage s;

    // Напечатать пустое хранилище без полей
    s.print();

    // Добавить скалярное поле int
    auto flag = s.add<int>("flag");

    // Векторное поле с 3 компонентами double
    auto vec = s.add<double>("vec", 3);

    // Рассматривается как пользовательский тип std::array<int, 3>,
    // то есть тип std::array<int, 3> не разворачивается
    auto arr = s.add<std::array<int, 3>>("arr");

    // Напечатать пустое хранилище, но есть имена полей
    s.print();

    // Увеличим размер до N элементов
    s.resize(9);

    // Типы инициализируются неизвестно чем
    s.print();

    {
        // Другой вариант, создать сразу хранилище на 100 элементов,
        // но без полей
        SoaStorage s2(100);
        s2.print();
    }

    // Заполнение
    for (size_t i = 0; i < s.size(); ++i) {
        // Скалярный тип int
        s(flag)[i] = i;

        // Аналогичное обращение через
        s.data(flag)[i] = i;

        // Векторный тип double
        s(vec[0])[i] = i + 0.1;
        s(vec[1])[i] = i + 0.2;
        s(vec[2])[i] = i + 0.3;

        // Пользовательский std::array<int, 3>
        s(arr)[i][0] = 100 + i;
        s(arr)[i][1] = 200 + i;
        s(arr)[i][2] = 300 + i;
    }

    // Пользовательские типы не выводятся, поскольку не известно,
    // как их интерпретировать
    s.print();

    // Выведем самостоятельно
    std::cout << "std::array<int, 3>:\n";
    for (size_t i = 0; i < s.size(); ++i) {
        std::cout << "    {" << s(arr)[i][0] << ", " << s(arr)[i][1] << " " << s(arr)[i][2] << "}\n";
    }
    std::cout << "\n";

    // Проверка наличия поля
    std::cout << "Contain int    'flag'? " << yes_no( s.contain<int>("flag") ) << " ( yes )" << "\n";
    std::cout << "Contain double 'flag'? " << yes_no( s.contain<double>("flag") ) << "  ( no  )" << "\n";

    // Восстановить переменную для доступа, если потеряли
    Storable<int> flag2 = s.storable<int>("flag");

    for (size_t i = 0; i < s.size(); ++i) {
        s(flag2)[i] = i * i;
    }
    // Массив "flag", естественно, изменился
    s.print();

    // Хранение сложных структур
    State z;

    // Изменяет хранилище s, добавляет в него новые поля
    z.add_types(s);

    // Также N элементов, но добавлены новые поля
    s.print();

    // Получить переменные для доступа к определенным полям,
    // нужно знать, как эти поля именует State
    auto rho = s.storable<double>("curr.rho");
    auto v   = s.storable<Vector3d>("half.v");

    // Лучше спросить у структуры напрямую
    v = z.half.v;

    for (size_t i = 0; i < s.size(); ++i) {
        s(z.curr.rho)[i] = i * i;
        s(z.next.v  )[i] = Vector3d{i + 0.1, i + 0.2, i + 0.3};

        s(z.half.P)  [i] = i * 1.0e5;
        s(z.next.rho)[i] = i * 1.1;
    }

    s.print();


    // Проверяем опцию копирования

    for (int i = 0; i < s.size(); ++i) {
        s(arr)[i] = {1 + i, 3 + i, 5 + i};
    }

    std::cout << "arr: \n";
    for (size_t i = 0; i < s.size(); ++i) {
        std::cout << "  " << i << ": " << s(arr)[i][0] << ", " << s(arr)[i][1] << ", " << s(arr)[i][2] << "\n";
    }

    s.copy_data(1, 6);

    s.print();

    std::cout << "arr: \n";
    for (size_t i = 0; i < s.size(); ++i) {
        std::cout << "  " << i << ": " << s(arr)[i][0] << ", " << s(arr)[i][1] << ", " << s(arr)[i][2] << "\n";
    }

    return 0;
}
