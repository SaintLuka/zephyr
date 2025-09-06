#include <iostream>
#include <iomanip>

#include <zephyr/mesh/storage.h>

using namespace zephyr::mesh;
using namespace zephyr::utils;
using zephyr::geom::Vector3d;

// Имитация вектора состояния
struct PState {
    Storable<double>   rho;  // Плотность
    Storable<Vector3d> v;    // Скорость
    Storable<double>   P;    // Давление
    Storable<double>   e;    // Энергия

    // Изменяет хранилище s, добавляет в него новые поля
    void add_types(Storage &s, const std::string &prefix) {
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

    // @brief Градиент вектора состояния
    PState d_dx, d_dy, d_dz;

    // Изменяет хранилище s, добавляет в него новые поля
    void add_types(Storage &s) {
        curr.add_types(s, "curr");
        half.add_types(s, "half");
        next.add_types(s, "next");
        d_dx.add_types(s, "d_dx");
        d_dy.add_types(s, "d_dy");
        d_dz.add_types(s, "d_dz");
    }
};

int main() {
    Storage s;

    // Добавить скалярное поле int
    auto flag = s.add<int>("flag");

    // Векторное поле с 3 компонентами double
    auto arr = s.add<double>("arr", 3);

    // Не эквивалентно предыдущему, Vector3d рассматривается как скалярный тип
    auto vec = s.add<Vector3d>("vec");

    // Напечатать пустое хранилище, но есть имена полей
    s.print(flag);
    s.print(arr);
    s.print(vec);

    // Увеличим размер до N элементов
    s.resize(9);

    // Типы инициализируются неизвестно чем
    s.print(flag);
    s.print(arr);
    s.print(vec);

    {
        // Создать сразу хранилище на 100 элементов, но без полей
        Storage s2(100);

        // Добавляется сразу массив на 100 элементов типа int
        s.add<int>("some");
    }

    // Заполнение
    for (size_t i = 0; i < s.size(); ++i) {
        // Скалярный тип int, возвращается ссылка
        s.get_val(flag, i) = i;

        // Векторный тип double, возвращается span

        // Можно записать в span из std::array или std::vector
        std::array<double, 3> some_arr = {i + 0.1, i + 0.2, i + 0.3};
        s.get_val(arr, i) = some_arr;

        // Можно обращаться по отдельным индексам
        s.get_val(arr, i)[2] = 0.0;

        // Класс Vector3d, возвращается ссылка Vector3d
        s.get_val(vec, i) = Vector3d{i + 0.5, i + 0.6, i + 0.7};
    }

    // Вывести можно только типы, у которых определен operator<<
    s.print(flag);
    s.print(arr);
    s.print(vec);

    // Проверяем опцию копирования
    s.copy_data(1, 6);

    std::cout << "\nAfter copy:\n";
    s.print(flag);
    s.print(arr);
    s.print(vec);

    // Проверка наличия поля
    std::cout << std::boolalpha << "\n";
    std::cout << "contain(\"flag\"):        " << s.contain("flag") << "  // = true\n";
    std::cout << "contain<short>(\"flag\"): " << s.contain<short>("flag") << " // = false\n";
    std::cout << "contain<float>(\"flag\"): " << s.contain<float>("flag") << "  // = true, scalar type and the same sizeof\n\n";

    // Восстановить переменную для доступа, если потеряли
    Storable<int> flag2 = s.find<int>("flag");

    for (size_t i = 0; i < s.size(); ++i) {
        s.get_val(flag2, i) = i * i;
    }
    // Массив "flag", естественно, изменился
    s.print(flag);

    // Хранение сложных структур
    State z;

    // Изменяет хранилище s, добавляет в него новые поля
    z.add_types(s);

    // Вывести названия новых полей
    std::cout << "\n";
    s.print_names();
    std::cout << "\n";

    // Получить переменные для доступа к определенным полям,
    // нужно знать, как эти поля именует State
    auto rho = s.find<double>("curr.rho");
    auto v   = s.find<Vector3d>("half.v");

    // Лучше спросить у структуры напрямую
    auto P = z.half.P;

    for (size_t i = 0; i < s.size(); ++i) {
        s.get_val(rho, i) = i * i;
        s.get_val(v  , i) = Vector3d{i + 0.1, i + 0.2, i + 0.3};
        s.get_val(P  , i) = i * 1.0e5;
        s.get_val(z.next.rho, i) = i * 1.1;
    }

    s.print(z.curr.rho);
    s.print(z.half.v);
    s.print(z.half.P);
    s.print(z.next.rho);
    return 0;
}