---
sidebar_position: 1
---

# Работа с сеткой

## Создание сетки

Для создания сетки используются сеточные генераторы. Наиболее простые генераторы `Rectangle` и `Cuboid` расположены в пространстве имён `zephyr::geom::generator`. Создание сетки:
```cpp
Rectangle gen(0.0, 1.0, 0.0, 0.5);  // [x_min, x_max, y_min, y_max]
gen.set_nx(200);
gen.set_boundaries({.left  =Boundary::ZOE, .right=Boundary::ZOE,
                    .bottom=Boundary::ZOE, .top  =Boundary::ZOE});

EuMesh mesh(gen);
```
Вызов `gen.set_nx(200)` устанавливает 200 ячеек по оси $x$, при этом число ячеек вдоль оси $y$ выбирается автоматически. Аналогично можно вызвать функцию `set_ny`, или же установить желаемое число ячеек по обеим осям:
```cpp
Rectangle gen(0.0, 1.0, 0.0, 0.5);
gen.set_sizes(200, 1);
```
В данном случае получится сетка с ячейками, вытянутыми вдоль оси $y$.

Генератор `Cuboid` работает аналогичным образом:
```cpp
Cuboid gen(0.0, 1.0, 0.0, 0.5, 0.0, 1.0);  // [x_min, x_max, y_min, y_max, z_min, z_max]
gen.set_nx(20);
gen.set_boundaries({.left  =Boundary::ZOE, .right=Boundary::ZOE,
                    .bottom=Boundary::ZOE, .top  =Boundary::ZOE,
                    .back  =Boundary::ZOE, .front=Boundary::ZOE});

EuMesh mesh(gen);
```

## Данные на сетке

Данные хранятся в ячейках сетки. Добавление новых полей данных выполняется командой `mesh.add<T>(name)`:
```cpp
Storable<double> u = mesh.add<double>("u");
```

На сетку можно добавить данные любого достаточно простого типа (тривиальный конструктор, простое копирование, отсутствие наследования). Переменная `u` в примере выше является "ключом", по которому можно получить данные ячейки.

Пройдем по всем ячейкам сетки и запишем значения по ключу `u`, доступ к данным выполняется через оператор `[]`:
```cpp
for (auto cell: mesh) {
  cell[u] = std::sin(cell.x());
}
```

Фактически, ключ `Storable<T>` это простой индекс:
```cpp
template<typename T>
struct Storable {
  int index = -1;
}
```
Данный индекс указывает на номер буфера данных в хранилище. Ключ `Storable<T>` можно копировать в произвольном количестве экземпляров. Поскольку он весит как `int`, его следует передавать в функции по значению, и по значению захватывать в лямбды.

## Цикл по ячейкам

Для обхода ячеек сетки предусмотрен range-based цикл. Также для каждой ячейки можно выполнить обход по граням. Для примера, найдем максимальное расстояние между парой смежных ячеек:
```cpp
double max_dist = 0.0;
for (auto cell: mesh) {
  Vector3d vc = cell.center();
  for (auto face: cell.faces()) {
    auto neib = face.neib();
    Vector3d vn = neib.center();
    double dist = (vc - vn).norm();
    max_dist = std::max(dist, max_dist);
  }
}
```

Через каждую грань можно получить соседнюю ячейку.

## Функции ячеек

|Фунция| Описание |
|-|-|
| `index() -> int` | Индекс ячейки на сетке |
| `level() -> int` | Уровень адаптации ячейки (с нуля) |
| `center() -> Vector3d` | Центр ячейки |
| `x(), y(), z() -> double` | Отдельные координаты центра ячейки |
| `volume() -> double` | Объем ячейки (в 3D) или площадь (в 2D) |
| `linear_size() -> double` | Линейные размеры ячейки |
| `incircle_diameter() -> double` | Диаметр вписанной окружности |
| `hx(), hy(), hz() -> double` | Размеры прямоугольной ячейки |
| `polygon() -> geom::Polygon` | Геометрия ячейки в виде двумерного полигона (не определено в 3D) |
| `polyhedron() -> geom::Polyhedron` | Геометрия ячейки в виде многогранника (не определено в 2D) |
| `simple_face(Side2D side) -> bool` | Простая грань по стороне? Пример: `cell.simple_face(Side2D::LEFT)`. |
| `simple_face(Side3D side) -> bool` | Простая грань по стороне? Пример: `cell.simple_face(Side3D::BACK)`. |
| `complex_face(Side2D side) -> bool` | Сложная грань по стороне? Пример: `cell.complex_face(Side2D::LEFT)`. |
| `complex_face(Side3D side) -> bool` | Сложная грань по стороне? Пример: `cell.complex_face(Side3D::BACK)`. |
| `face(Side2D side) -> EuFace` | Получить грань ячейки: `cell.face(Side2D::LEFT)`, подгрань `Side2D::LEFT[1]`. |
| `face(Side3D side) -> EuFace` | Получить грань ячейки: `cell.face(Side3D::BACK)`, подгрань `Side3D::BACK[3]`. |
| `faces(Direction dir = ANY) -> EuFaces` | Получить итератор по граням. |

Итератор по граням позволяет ограничить обход только одним направлением
```cpp
for (auto face: cell.faces(Direction::X)) {
  // ..
}
```
В данном случае обход осуществляется только по граням с нормалью вдоль оси x.


## Функции граней

| Функция | Описание |
|-|-|
| `flag() -> Boundary` | Флаг граничных условий. |
| `is_bondary() -> bool` | Грань на границе? Флаг `Boundary::PERIODIC` не считается граничным условием. |
| `normal() -> Vector3d` | Внешняя нормаль по отношению к исходной ячейке. |
| `center() -> Vector3d` | Центр грани. |
| `x(), y(), z() -> Vector3d`| Отдельные координаты центра грани. |
| `area() -> double` | Площадь грани (в 3D), длина грани (в 2D). |
| `area_n() -> Vector3d` | Площадь грани, умноженная на внешнюю нормаль. |
| `symm_point(Vector3d p) -> Vector3d` | Получить точку, симметричную относительно грани. |
| `neib() -> EuCell` | Соседняя ячейка через грань |

Функция `neib()` всегда актуальна. Если сосед отсутствует, тогда создается `EuCell`, которая ссылается на исходную ячейку.

Некоторые функции позволяют сразу получить данные соседней ячейки без создания дополнительного экземпляра `EuCell` для соседней ячейки.
```cpp
auto neib = face.neib();
Vector3d nc1 = neib.center();
Vector3d nc2 = face.neib_center();
```

Такие функции будут немного быстрее.

| | |
|-|-|
| `neib(Storable<T> key) -> const T&` | Получить данные из соседней ячейки. |
| `neib_flag() -> int` | Флаг адаптации соседней ячейки. |
| `neib_center() -> Vector3d` | Центр соседней ячейки. |
| `neib_volume() -> double` | Объем соседней ячейки. |





