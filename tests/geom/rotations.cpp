#include <zephyr/mesh/amr/rotations.h>

using namespace zephyr::mesh::amr;

template<typename Matrix>
void print_matrix(const Matrix& m) {
    std::cout << "{";
    for (int i = 0; i < m.rows(); ++i) {
        std::cout << "{";
        for (int j = 0; j < m.cols(); ++j) {
            std::cout << std::format("{:2d}", m(i, j));
            if (j + 1 < m.cols()) std::cout << ",";
        }
        std::cout << "}";
        if (i + 1 < m.rows()) std::cout << ", ";
    }
    std::cout << "}";
}

// Вывод матриц в виде constexpr-инициализации
inline void print_rotation_matrices() {
    auto m2 = generate_rotation_matrices_2d();
    auto m3 = generate_rotation_matrices_3d();

    std::cout << "static constexpr int rotation_matrices_2d[8][2][2] = {\n";
    for (const auto& m : m2) {
        std::cout << "  ";
        print_matrix(m);
        std::cout << ",\n";
    }
    std::cout << "};\n\n";

    std::cout << "static constexpr int rotation_matrices_3d[48][3][3] = {\n";
    for (const auto& m : m3) {
        std::cout << "  ";
        print_matrix(m);
        std::cout << ",\n";
    }
    std::cout << "};\n\n";
}

template<size_t N>
void print_perm(const std::array<int, N>& p) {
    std::cout << "{";
    for (size_t i = 0; i < N; ++i) {
        std::cout << p[i];
        if (i + 1 < N) std::cout << ", ";
    }
    std::cout << "}";
}

// Вывод перестановок в виде constexpr-инициализации
inline void print_permutations() {
    auto p2 = generate_permutations_2d();
    auto p3 = generate_permutations_3d();

    std::cout << "/// @brief Полная группа перестановок квадрата, всего 8 элементов.\n";
    std::cout << "/// Элементы записаны в форме перестановок индексов вершин квадрата\n";
    std::cout << "static constexpr int permutations_2d[8][4] = {\n";
    for (const auto& p : p2) {
        std::cout << "    ";
        print_perm(p);
        std::cout << ",\n";
    }
    std::cout << "};\n\n";

    std::cout << "/// @brief Полная группа перестановок куба, всего 48 элементов\n";
    std::cout << "/// Элементы записаны в форме перестановок индексов вершин куба.\n";
    std::cout << "static constexpr int permutations_3d[48][8] = {\n";
    for (const auto& p : p3) {
        std::cout << "    ";
        print_perm(p);
        std::cout << ",\n";
    }
    std::cout << "};\n\n";
}

inline void subface_by_side_child_gen() {
    int gen_subface_by_side_child_2D[4][4] = {
        {Side2D::L[0], -1, Side2D::L[1], -1},
        {-1, Side2D::R[0], -1, Side2D::R[1]},
        {Side2D::B[0], Side2D::B[1], -1, -1},
        {-1, -1, Side2D::T[0], Side2D::T[1]}
    };
    // Проходим по subfaces
    for (Side2D s = 0; s < Side2D::n_subfaces(); ++s.val) {
        if (gen_subface_by_side_child_2D[s.simple()][indexing::child(s)] != s) {

            std::cout << s.simple() << " " << indexing::child(s) << " " << s << "\n";

            throw std::runtime_error("No way 2D");
        }
    }

    std::cout << "/// @brief Получить subface по стороне и индексу дочерней ячейки [side][child]\n";
    std::cout << "static constexpr int subface_by_side_child_2D[4][4] = {\n";
    for (Side2D side: Side2D::items()) {
        std::cout << "    {";
        for (int i = 0; i < 3; ++i) {
            std::cout << gen_subface_by_side_child_2D[side][i] << ", ";
        }
        std::cout << gen_subface_by_side_child_2D[side][3] << "},\n";
    }
    std::cout << "};\n\n";

    int gen_subface_by_side_child_3D[6][8] = {
        {Side3D::L[0], -1, Side3D::L[1], -1, Side3D::L[2], -1, Side3D::L[3], -1},
        {-1, Side3D::R[0], -1, Side3D::R[1], -1, Side3D::R[2], -1, Side3D::R[3]},
        {Side3D::B[0], Side3D::B[1], -1, -1, Side3D::B[2], Side3D::B[3], -1, -1},
        {-1, -1, Side3D::T[0], Side3D::T[1], -1, -1, Side3D::T[2], Side3D::T[3]},
        {Side3D::Z[0], Side3D::Z[1], Side3D::Z[2], Side3D::Z[3], -1, -1, -1, -1},
        {-1, -1, -1, -1, Side3D::F[0], Side3D::F[1], Side3D::F[2], Side3D::F[3]},
    };
    // Проходим по subfaces
    for (Side3D s = 0; s < Side3D::n_subfaces(); ++s.val) {
        if (gen_subface_by_side_child_3D[s.simple()][indexing::child(s)] != s) {
            throw std::runtime_error("No way 3D");
        }
    }
    std::cout << "/// @brief Получить subface по стороне и индексу дочерней ячейки [side][child]\n";
    std::cout << "static constexpr int subface_by_side_child_3D[6][8] = {\n";
    for (Side3D side: Side3D::items()) {
        std::cout << "    {";
        for (int i = 0; i < 7; ++i) {
            std::cout << gen_subface_by_side_child_3D[side][i] << ", ";
        }
        std::cout << gen_subface_by_side_child_3D[side][7] << "},\n";
    }
    std::cout << "};\n\n";
}

inline void neib_child_gen_2D() {
    static const auto group = generate_permutations_2d();
    int gen_neib_child_by_symm_subface_2D[group.size()][Side2D::n_subfaces()];

    std::array<Vector3d, 4> UQ = {
        Vector3d{-0.5, -0.5, 0.0},
        Vector3d{+0.5, -0.5, 0.0},
        Vector3d{-0.5, +0.5, 0.0},
        Vector3d{+0.5, +0.5, 0.0}
    };

    SqQuad quad(UQ[0], UQ[1], UQ[2], UQ[3]);

    for (int r = 0; r < group.size(); ++r) {
        auto RQ = permute_2d(UQ, r);

        std::array<SqQuad, Side2D::count()> neibs;
        for (Side2D side: Side2D::items()) {
            neibs[side] = SqQuad(RQ[0], RQ[1], RQ[2], RQ[3]);
        }

        for (int i = 0; i < 9; ++i) {
            neibs[Side2D::L][i].x() -= 1.0;
            neibs[Side2D::R][i].x() += 1.0;
            neibs[Side2D::B][i].y() -= 1.0;
            neibs[Side2D::T][i].y() += 1.0;
        }

        std::array<std::array<Vector3d, indexing::CpC(2)>, Side2D::count()> child_c;
        for (Side2D side: Side2D::items()) {
            child_c[side] = {
                neibs[side](-0.5, -0.5),
                neibs[side](+0.5, -0.5),
                neibs[side](-0.5, +0.5),
                neibs[side](+0.5, +0.5)
            };
        }

        for (Side2D side = 0; side < Side2D::n_subfaces(); ++side.val) {
            Vector3d fc = 0.5 * (quad[indexing::amr::cf(side)[0]] + quad[indexing::amr::cf(side)[1]]);

            double min_dist = 100.0;
            int child_idx = -1;
            for (int i = 0; i < indexing::CpC(2); ++i) {
                double dist = (child_c[side.simple()][i] - fc).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    child_idx = i;
                }
            }
            gen_neib_child_by_symm_subface_2D[r][side] = child_idx;
        }
    }

    std::cout << "/// @brief Индекс дочерней ячейки соседа [symm][subface]\n";
    std::cout << "static constexpr int neib_child_by_symm_subface_2D[8][8] = {\n";
    for (int i = 0; i < group.size(); ++i) {
        std::cout << "    {";
        for (int j = 0; j < Side2D::n_subfaces() - 1; ++j) {
            std::cout << gen_neib_child_by_symm_subface_2D[i][j] << ", ";
        }
        std::cout << gen_neib_child_by_symm_subface_2D[i][Side2D::n_subfaces() - 1] << "},\n";
    }
    std::cout << "};\n\n";
}

inline void neib_child_gen_3D() {
    static const auto group = generate_permutations_3d();

    int gen_neib_child_by_symm_subface_3D[group.size()][Side3D::n_subfaces()];

    std::array<Vector3d, 8> UQ = {
        Vector3d{-0.5, -0.5, -0.5},
        Vector3d{+0.5, -0.5, -0.5},
        Vector3d{-0.5, +0.5, -0.5},
        Vector3d{+0.5, +0.5, -0.5},
        Vector3d{-0.5, -0.5, +0.5},
        Vector3d{+0.5, -0.5, +0.5},
        Vector3d{-0.5, +0.5, +0.5},
        Vector3d{+0.5, +0.5, +0.5}
    };

    SqCube quad(UQ[0], UQ[1], UQ[2], UQ[3], UQ[4], UQ[5], UQ[6], UQ[7]);

    for (int r = 0; r < group.size(); ++r) {
        auto RQ = permute_3d(UQ, r);

        std::array<SqCube, Side3D::count()> neibs;
        for (Side3D side: Side3D::items()) {
            neibs[side] = SqCube(RQ[0], RQ[1], RQ[2], RQ[3], RQ[4], RQ[5], RQ[6], RQ[7]);
        }

        for (int i = 0; i < 27; ++i) {
            neibs[Side3D::L][i].x() -= 1.0;
            neibs[Side3D::R][i].x() += 1.0;
            neibs[Side3D::B][i].y() -= 1.0;
            neibs[Side3D::T][i].y() += 1.0;
            neibs[Side3D::Z][i].z() -= 1.0;
            neibs[Side3D::F][i].z() += 1.0;
        }

        std::array<std::array<Vector3d, indexing::CpC(3)>, Side3D::count()> child_c;
        for (Side3D side: Side3D::items()) {
            child_c[side] = {
                neibs[side](-0.5, -0.5, -0.5),
                neibs[side](+0.5, -0.5, -0.5),
                neibs[side](-0.5, +0.5, -0.5),
                neibs[side](+0.5, +0.5, -0.5),
                neibs[side](-0.5, -0.5, +0.5),
                neibs[side](+0.5, -0.5, +0.5),
                neibs[side](-0.5, +0.5, +0.5),
                neibs[side](+0.5, +0.5, +0.5)
            };
        }

        for (Side3D side = 0; side < Side3D::n_subfaces(); ++side.val) {
            Vector3d fc = 0.25 * (
                quad[indexing::amr::cf(side)[0]] + quad[indexing::amr::cf(side)[1]] +
                quad[indexing::amr::cf(side)[2]] + quad[indexing::amr::cf(side)[3]]);

            double min_dist = 100.0;
            int child_idx = -1;
            for (int i = 0; i < indexing::CpC(3); ++i) {
                double dist = (child_c[side.simple()][i] - fc).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    child_idx = i;
                }
            }
            gen_neib_child_by_symm_subface_3D[r][side] = child_idx;
        }
    }

    std::cout << "/// @brief Индекс дочерней ячейки соседа [symm][subface]\n";
    std::cout << "static constexpr int neib_child_by_symm_subface_3D[48][24] = {\n";
    for (int i = 0; i < group.size(); ++i) {
        std::cout << "    {";
        for (int j = 0; j < Side3D::n_subfaces() - 1; ++j) {
            std::cout << gen_neib_child_by_symm_subface_3D[i][j] << ", ";
        }
        std::cout << gen_neib_child_by_symm_subface_3D[i][Side3D::n_subfaces() - 1] << "},\n";
    }
    std::cout << "};\n\n";
}

Eigen::Matrix2d random_rotation_2d(double max_angle) {
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution angle_dist(-max_angle, max_angle);

    // Случайный угол
    double phi = angle_dist(gen);

    Matrix2d R;
    R << std::cos(phi), std::sin(phi), -std::sin(phi), std::cos(phi);
    return R;
}

Matrix3d random_rotation_3d(double max_angle) {
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution angle_dist(-max_angle, max_angle);
    std::uniform_real_distribution axis_dist(-1.0, 1.0);

    // Случайная ось
    Vector3d axis(axis_dist(gen), axis_dist(gen), axis_dist(gen));
    axis.normalize();

    // Случайный угол
    double phi = angle_dist(gen);

    return Eigen::AngleAxisd(phi, axis).matrix();
}

void test_2d() {
    // Взяли квадрат и случайно его повернули
    Quad quad1(Vector3d{-1, -1, 0.0},
               Vector3d{+1, -1, 0.0},
               Vector3d{-1, +1, 0.0},
               Vector3d{+1, +1, 0.0});

    // Сильный поворот
    Matrix2d R1 = random_rotation_2d(3.15);
    for (int j = 0; j < 4; ++j) {
        quad1[j].head<2>() = R1 * quad1[j].head<2>();
    }

    // Оси квадрата
    Vector2d ex = (quad1.get(1.0, 0.0) - quad1.get(-1.0, 0.0)).head<2>().normalized();
    Vector2d ey = (quad1.get(0.0, 1.0) - quad1.get(0.0, -1.0)).head<2>().normalized();

    auto group = generate_rotation_matrices_2d();

    for (int r = 0; r < group.size(); ++r) {
        // Создали второй квадрат, слабо повернули на случайный угол
        Quad quad2(quad1);
        Matrix2d R;
        for (int j = 0; j < 4; ++j) {
            R.data()[j] = group[r].data()[j];
        }
        R = R * random_rotation_2d(0.25);
        for (int j = 0; j < 4; ++j) {
            quad2[j].head<2>() = R * quad2[j].head<2>();
        }

        // Оси куба
        Vector2d ex2 = (quad2.get(1.0, 0.0) - quad2.get(-1.0, 0.0)).head<2>().normalized();
        Vector2d ey2 = (quad2.get(0.0, 1.0) - quad2.get(0.0, -1.0)).head<2>().normalized();

        int r2 = find_rotation(axes2d_t{ex, ey}, axes2d_t{ex2, ey2});
        if (r != r2) {
            throw std::runtime_error("Can't define rotation 2D");
        }
    }
    std::cout << "Test 2D: OK\n";
}

void test_3d() {
    // Взяли куб и случайно его повернули
    Cube cube1( Vector3d{-1, -1, -1},
                Vector3d{+1, -1, -1},
                Vector3d{-1, +1, -1},
                Vector3d{+1, +1, -1},
                Vector3d{-1, -1, +1},
                Vector3d{+1, -1, +1},
                Vector3d{-1, +1, +1},
                Vector3d{+1, +1, +1});

    // Сильный поворот
    Matrix3d R1 = random_rotation_3d(3.15);
    for (int j = 0; j < 8; ++j) {
        cube1[j] = R1 * cube1[j];
    }

    // Оси куба
    Vector3d ex = (cube1.get(1.0, 0.0, 0.0) - cube1.get(-1.0, 0.0, 0.0)).normalized();
    Vector3d ey = (cube1.get(0.0, 1.0, 0.0) - cube1.get(0.0, -1.0, 0.0)).normalized();
    Vector3d ez = (cube1.get(0.0, 0.0, 1.0) - cube1.get(0.0, 0.0, -1.0)).normalized();

    auto group = generate_rotation_matrices_3d();

    for (int r = 0; r < group.size(); ++r) {
        // Создали второй куб, слабо повернули на случайный угол
        Cube cube2(cube1);
        Matrix3d R;
        for (int j = 0; j < 9; ++j) {
            R.data()[j] = group[r].data()[j];
        }

        R = R * random_rotation_3d(0.25);
        for (int j = 0; j < 8; ++j) {
            cube2[j] = R * cube2[j];
        }

        // Оси куба
        Vector3d ex2 = (cube2.get(1.0, 0.0, 0.0) - cube2.get(-1.0, 0.0, 0.0)).normalized();
        Vector3d ey2 = (cube2.get(0.0, 1.0, 0.0) - cube2.get(0.0, -1.0, 0.0)).normalized();
        Vector3d ez2 = (cube2.get(0.0, 0.0, 1.0) - cube2.get(0.0, 0.0, -1.0)).normalized();

        int r2 = find_rotation(axes3d_t{ex, ey, ez}, axes3d_t{ex2, ey2, ez2});
        if (r != r2) {
            throw std::runtime_error("Can't define rotation 3D");
        }
    }
    std::cout << "Test 3D: OK\n";
}

int main() {
    //print_rotation_matrices();
    //print_permutations();
    subface_by_side_child_gen();
    neib_child_gen_2D();
    neib_child_gen_3D();

    test_2d();
    test_3d();

    return 0;
}