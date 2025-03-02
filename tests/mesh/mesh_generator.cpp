#include <iostream>
#include <vector>

#include <zephyr/mesh/mesh.h>
#include <zephyr/io/vtu_file.h>

#include <zephyr/geom/box.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/cubic.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/block_structured.h>


using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::mesh::generator;

// Данные ячеек
struct _U_ {
    int uid;
    double val;
};

_U_ U;

// Данные узлов
struct _V_ {
    int uid;
    double val;
};

_V_ V;

double get_uid_amr(AmrStorage::Item& cell) {
    return cell(U).uid;
}

double get_val_amr(AmrStorage::Item& cell) {
    return cell(U).val;
}

double get_uid_cell(CellStorage::Item& cell) {
    return cell(U).uid;
}

double get_val_cell(CellStorage::Item& cell) {
    return cell(U).val;
}

double get_uid_node(NodeStorage::Item& cell) {
    return cell(V).uid;
}

double get_val_node(NodeStorage::Item& cell) {
    return cell(V).val;
}

struct Test {
    enum TestType {
        Rectangle, RectangleHex, Cuboid,
        Sector1, Sector2, Sector3, Disk,
        BlockStruct1, BlockStruct2, BlockStruct3
    };

    std::string name;  ///< Название теста
    std::string file;  ///< Имя выходного файла
    std::string desc;  ///< Описание теста
    
    Generator::Ptr generator;  ///< Сеточный генератор
    
    Test(TestType test);

    EuMesh gen_eu() const;

    LaMesh gen_la() const;
};

int main() {
    // Переменные для записи
    Variables vars;
    vars += { "uid", get_uid_amr };
    vars += { "val", get_val_amr };
    vars += { "uid", get_uid_cell };
    vars += { "val", get_val_cell };
    vars += { "uid", get_uid_node };
    vars += { "val", get_val_node };

    // Список тестов
    std::vector<Test> tests_list = {
            Test::Rectangle,
            Test::RectangleHex,
            Test::Cuboid,
            Test::Sector1,
            Test::Sector2,
            Test::Sector3,
            Test::Disk,
            Test::BlockStruct1,
            Test::BlockStruct2,
            Test::BlockStruct3
    };

    for (auto& test: tests_list) {
        std::cout << "Тест '" << test.name << "'.\n";
        std::cout << "\t" << test.desc << "\n";
        std::cout << "\tВыходной файл: " << test.file << "\n\n";

        VtuFile file("output/" + test.file, vars);
        file.hex_only = false;
        file.unique_nodes = true;

        auto cells = test.gen_eu();
        //auto cells = test.gen_la();
        file.save(cells);
    }
}

Generator::Ptr create_rectangle() {
    auto rect = Rectangle::create(0.0, 2.0, 0.0, 1.0, false);
    rect->set_nx(1000);
    return rect;
}

Generator::Ptr create_rectangle_hex() {
    auto rect = Rectangle::create(0.0, 2.0, 0.0, 1.0, true);
    rect->set_nx(100);
    return rect;
}

Generator::Ptr create_cuboid() {
    auto cuboid = Cuboid::create(0.0, 2.0, 0.0, 1.0, 0.0, 0.5);
    cuboid->set_nx(40);
    return cuboid;
}

Generator::Ptr create_sector1() {
    auto sector = Sector::create(1.0, 0.3, M_PI / 5.0, true);
    sector->set_n_phi(17);
    return sector;
}

Generator::Ptr create_sector2() {
    auto sector = Sector::create(1.0, 0.3, 0.8 * M_PI, false);
    sector->set_n_phi(28);
    return sector;
}

Generator::Ptr create_sector3() {
    auto sector = Sector::create(1.0, 0.3, 1.4 * M_PI, false);
    sector->set_n_phi(36);
    return sector;
}

Generator::Ptr create_disk() {
    auto sector = Sector::create(1.0, 0.3, 2.0 * M_PI, false);
    sector->set_n_phi(80);
    return sector;
}

Generator::Ptr create_block_structured1() {
    double R = 1.0;

    // Задаем базисные вершины для струтурированных блоков
    BaseVertex::Ptr v1 = BaseVertex::create(0.0, 0.0, true);
    BaseVertex::Ptr v2 = BaseVertex::create(R / 2, 0.0, false);
    BaseVertex::Ptr v3 = BaseVertex::create(R, 0.0, true);
    BaseVertex::Ptr v4 = BaseVertex::create(0.0, R / 2, false);

    BaseVertex::Ptr v5 = BaseVertex::create(R / 2, R / 2, false);
    BaseVertex::Ptr v6 = BaseVertex::create(R / sqrt(2.0), R / sqrt(2.0), false);
    BaseVertex::Ptr v7 = BaseVertex::create(0.0, R, true);

    // Дополнительные вершины для задания кубического сплайна
    Vector3d v_s1 = {R/3, 0.0, 0.0};
    Vector3d v_s2 = {R/2, R/9, 0.0};

    // Ограничивающие кривые области
    Curve::Ptr circle = Circle::create(v3, v6, v7);           // Окружность
    Curve::Ptr left = Plane::create(v1, v7);                   // Прямая слева
    Curve::Ptr bottom = Cubic::create(*v1, v_s1, v_s2, *v3);  // Сплайн на нижней границе

    // Генератор сетки
    auto gen = BlockStructured::create(3);
    BlockStructured& blocks = *gen;

    blocks[0] = {v1, v2, v4, v5};
    blocks[0].set_boundary(v1, v4, left);
    blocks[0].set_boundary(v1, v2, bottom);

    blocks[1] = {v2, v3, v5, v6};
    blocks[1].set_boundary(v2, v3, bottom);
    blocks[1].set_boundary(v3, v6, circle);

    blocks[2] = {v4, v5, v6, v7};
    blocks[2].set_boundary(v4, v7, left);
    blocks[2].set_boundary(v6, v7, circle);

    // Необходимо связать блоки
    blocks.link();

    // Характерный размер ячейки
    double h = 0.02;

    // Число ячеек по грани в блоках
    size_t N = size_t(0.5 * R / h);

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    blocks[0].set_size(v1, v2, N);
    blocks[0].set_size(v1, v4, N);
    blocks[1].set_size(v2, v3, N);

    // Точность сглаживания
    blocks.set_accuracy(1.0e-5);

    // Возвращаем генератор
    return gen;
}

Generator::Ptr create_block_structured2() {
    double L = 40.0;
    double H = 10.0;
    double r = 0.8;
    double R = 3.5;
    double a = r / std::sqrt(2.0);

    // Задаем базисные вершины для струтурированных блоков
    BaseVertex::Ptr v1 = BaseVertex::create(-L, -H, true);
    BaseVertex::Ptr v2 = BaseVertex::create(-R, -H, false);
    BaseVertex::Ptr v3 = BaseVertex::create(+R, -H, false);
    BaseVertex::Ptr v4 = BaseVertex::create(+L, -H, true);

    BaseVertex::Ptr v5 = BaseVertex::create(-L, -R, false);
    BaseVertex::Ptr v6 = BaseVertex::create(-R, -R, false);
    BaseVertex::Ptr v7 = BaseVertex::create(+R, -R, false);
    BaseVertex::Ptr v8 = BaseVertex::create(+L, -R, false);

    BaseVertex::Ptr v9  = BaseVertex::create(-a, -a, false);
    BaseVertex::Ptr v10 = BaseVertex::create(+a, -a, false);
    BaseVertex::Ptr v11 = BaseVertex::create(-a, +a, false);
    BaseVertex::Ptr v12 = BaseVertex::create(+a, +a, false);

    BaseVertex::Ptr v13 = BaseVertex::create(-L, +R, false);
    BaseVertex::Ptr v14 = BaseVertex::create(-R, +R, false);
    BaseVertex::Ptr v15 = BaseVertex::create(+R, +R, false);
    BaseVertex::Ptr v16 = BaseVertex::create(+L, +R, false);

    BaseVertex::Ptr v17 = BaseVertex::create(-L, +H, true);
    BaseVertex::Ptr v18 = BaseVertex::create(-1.5*R, +1.3*H, false);
    BaseVertex::Ptr v19 = BaseVertex::create(+1.5*R, +0.7*H, false);
    BaseVertex::Ptr v20 = BaseVertex::create(+L, +H, true);

    // Ограничивающие кривые области
    Curve::Ptr circle = Circle::create(v9, v10, v11);
    Curve::Ptr left   = Plane::create(v1, v17);
    Curve::Ptr right  = Plane::create(v4, v20);
    Curve::Ptr bottom = Plane::create(v1, v4);
    Curve::Ptr top    = Cubic::create(v17, v18, v19, v20);

    // Генератор сетки
    auto gen = BlockStructured::create(12);
    BlockStructured& blocks = *gen;

    blocks[0] = {v1, v2, v5, v6};
    blocks[0].set_boundary(v1, v5, left);
    blocks[0].set_boundary(v1, v2, bottom);

    blocks[1] = {v2, v3, v7, v6};
    blocks[1].set_boundary(v2, v3, bottom);

    blocks[2] = {v3, v4, v8, v7};
    blocks[2].set_boundary(v3, v4, bottom);
    blocks[2].set_boundary(v4, v8, right);

    blocks[3] = {v5, v6, v14, v13};
    blocks[3].set_boundary(v5, v13, left);

    blocks[4] = {v6, v9, v11, v14};
    blocks[4].set_boundary(v9, v11, circle);

    blocks[5] = {v6, v7, v10, v9};
    blocks[5].set_boundary(v9, v10, circle);

    blocks[6] = {v10, v7, v15, v12};
    blocks[6].set_boundary(v10, v12, circle);

    blocks[7] = {v11, v12, v15, v14};
    blocks[7].set_boundary(v11, v12, circle);

    blocks[8] = {v7, v8, v16, v15};
    blocks[8].set_boundary(v8, v16, right);

    blocks[9] = {v13, v14, v18, v17};
    blocks[9].set_boundary(v13, v17, left);
    blocks[9].set_boundary(v17, v18, top);

    blocks[10] = {v14, v15, v19, v18};
    blocks[10].set_boundary(v18, v19, top);

    blocks[11] = {v15, v16, v20, v19};
    blocks[11].set_boundary(v16, v20, right);
    blocks[11].set_boundary(v19, v20, top);

    // Необходимо связать блоки
    blocks.link();

    // Характерный размер ячейки
    double DX = 0.5;

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    blocks[0].set_size(v1, v2, size_t((L - R) / DX));
    blocks[1].set_size(v2, v3, size_t(2 * R / DX));
    blocks[2].set_size(v3, v4, size_t((L - R) / DX));
    blocks[10].set_size(v18, v19, size_t(2 * R / DX));

    blocks[0].set_size(v1, v5, size_t((H - R) / DX));
    blocks[3].set_size(v5, v13, size_t(2 * R / DX));
    blocks[9].set_size(v13, v17, size_t((H - R) / DX));
    blocks[8].set_size(v8, v16, size_t(2 * R / DX));

    blocks[4].set_size(v6, v9, size_t(1.5*(R * std::sqrt(2.0) - r) / DX));

    // Точность сглаживания (необязательно)
    blocks.set_accuracy(1.0e-5);

    // Возвращаем генератор
    return gen;
}

Generator::Ptr create_block_structured3() {
    double R = 8.0;
    double r = 0.5;
    double h = 1.0;

    // Задаем базисные вершины для струтурированных блоков
    BaseVertex::Ptr v1 = BaseVertex::create(0.0, 0.0, true);
    BaseVertex::Ptr v2 = BaseVertex::create(h - r, 0.0, true);
    BaseVertex::Ptr v3 = BaseVertex::create(h + r, 0.0, true);
    BaseVertex::Ptr v4 = BaseVertex::create(h + R, 0.0, false);
    BaseVertex::Ptr v5 = BaseVertex::create(0.0, r, false);
    BaseVertex::Ptr v6 = BaseVertex::create(h, r, false);
    BaseVertex::Ptr v7 = BaseVertex::create(0.0, R, true);
    BaseVertex::Ptr v8 = BaseVertex::create(h, R, false);

    // Ограничивающие кривые области
    Curve::Ptr circ = Circle::create(r, {h, 0.0, 0.0});
    Curve::Ptr CIRC = Circle::create(v7, v8, v4);
    Curve::Ptr left = Plane::create(v1, v7);
    Curve::Ptr bottom = Plane::create(v1, v4);

    // Генератор сетки
    auto gen = BlockStructured::create(3);
    BlockStructured& blocks = *gen;

    blocks[0] = {v1, v2, v6, v5};
    blocks[0].set_boundary(v1, v2, bottom);
    blocks[0].set_boundary(v1, v5, left);
    blocks[0].set_boundary(v2, v6, circ);

    blocks[1] = {v5, v6, v8, v7};
    blocks[1].set_boundary(v5, v7, left);
    blocks[1].set_boundary(v7, v8, CIRC);

    blocks[2] = {v3, v4, v8, v6};
    blocks[2].set_boundary(v6, v3, circ);
    blocks[2].set_boundary(v4, v8, CIRC);
    blocks[2].set_boundary(v3, v4, bottom);

    // Необходимо связать блоки
    blocks.link();

    // Характерный размер ячейки
    double DX = 0.1;

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    blocks[0].set_size(v1, v2, 1.5*(h - r) / DX);
    blocks[0].set_size(v1, v5, r / DX);

    blocks[2].set_size(v3, v6, 4 * r / DX);
    blocks[2].set_size(v3, v4, 0.95*log((R+h)/r)/log(1 + M_PI*DX/2*r));

    // Точность сглаживания
    blocks.set_accuracy(1.0e-5);

    // Возвращаем генератор
    return gen;
}

Test::Test(TestType test) {
    switch (test) {
        case Test::Rectangle:
            name = "Rectangle";
            file = "rectangle.vtu";
            desc = "Декартова сетка в прямоуольнике.";
            generator = create_rectangle();
            break;

        case Test::RectangleHex:
            name = "RectangleHex";
            file = "rectangle_hex.vtu";
            desc = "Сетка из шестиугольников в прямоугольнике.";
            generator = create_rectangle_hex();
            break;

        case Test::Cuboid:
            name = "Cuboid";
            file = "cuboid.vtu";
            desc = "Декартова сетка в прямоугольном параллелепипеде.";
            generator = create_cuboid();
            break;

        case Test::Sector1:
            name = "Sector 1";
            file = "sector_1.vtu";
            desc = "Сетка в секторе с выколотым центром.";
            generator = create_sector1();
            break;

        case Test::Sector2:
            name = "Sector 2";
            file = "sector_2.vtu";
            desc = "Сетка в секторе (угол раствора < п)";
            generator = create_sector2();
            break;

        case Test::Sector3:
            name = "Sector 3";
            file = "sector_3.vtu";
            desc = "Сетка в секторе (угол раствора > п)";
            generator = create_sector3();
            break;

        case Test::Disk:
            name = "Disk";
            file = "disk.vtu";
            desc = "Блочно-структурированная сетка в круге.";
            generator = create_disk();
            break;

        case Test::BlockStruct1:
            name = "Block Structured 1";
            file = "block_structured_1.vtu";
            desc = "Блочно-структурированная сетка №1.";
            generator = create_block_structured1();
            break;

        case Test::BlockStruct2:
            name = "Block Structured 2";
            file = "block_structured_2.vtu";
            desc = "Блочно-структурированная сетка №2.";
            generator = create_block_structured2();
            break;

        case Test::BlockStruct3:
            name = "Block Structured 3";
            file = "block_structured_3.vtu";
            desc = "Блочно-структурированная сетка №3.";
            generator = create_block_structured3();
            break;

        default:
            throw std::runtime_error("Unknown test");
    }
}

void fill(AmrStorage& cells) {
}

EuMesh Test::gen_eu() const {
    EuMesh mesh(generator, U);
    
    Box box = mesh.bbox();
    double L = box.diameter();

    for (auto& cell: mesh) {
        Vector3d r = cell.center();
        cell(U).uid = cell - mesh.begin();
        cell(U).val = std::cos(1.0 / (r.squaredNorm() / (L * L) + 0.03));
    }

    return mesh;
}

LaMesh Test::gen_la() const {
    LaMesh mesh(generator.get(), U, V);

    Box box = mesh.bbox();
    double L = box.diameter();
    
    for (auto& cell: mesh.locals()) {
        Vector3d r = cell.center;
        cell(U).uid = cell.index;
        cell(U).val = std::cos(1.0 / (r.squaredNorm() / (L * L) + 0.03));
    }

    for (auto& node: mesh.nodes()) {
        Vector3d r = node.coords;
        node(V).uid = node.index;
        node(V).val = std::sin(1.0 / (r.squaredNorm() / (L * L) + 0.03));
    }

    return mesh;
}
