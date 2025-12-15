#include <zephyr/geom/vector.h>
#include <zephyr/geom/geom.h>
#include <zephyr/mesh/side.h>

using namespace zephyr;
using namespace zephyr::geom;
using namespace zephyr::mesh;

// Полная группа перестановок квадрата, элементы записаны в форме перестановок
// индексов вершин квадрата.
static constexpr std::array<std::array<int, 4>, 8> quad_symmetries = {std::array<int, 4>
    {0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 0, 1}, {3, 0, 1, 2},
    {2, 1, 0, 3}, {3, 2, 1, 0}, {0, 3, 2, 1}, {1, 0, 3, 2},
};

// Полная группа перестановок куба.
// Всего 48 элементов, я их сгенерировал каким-то скриптом в python,
// элементы записаны в форме перестановок индексов вершин куба.
static constexpr std::array<std::array<int, 8>, 48> cube_symmetries = { std::array<int, 8>
    {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 4, 5, 2, 3, 6, 7}, {0, 2, 1, 3, 4, 6, 5, 7},
    {0, 2, 4, 6, 1, 3, 5, 7}, {0, 4, 1, 5, 2, 6, 3, 7}, {0, 4, 2, 6, 1, 5, 3, 7},
    {1, 0, 3, 2, 5, 4, 7, 6}, {1, 0, 5, 4, 3, 2, 7, 6}, {1, 3, 0, 2, 5, 7, 4, 6},
    {1, 3, 5, 7, 0, 2, 4, 6}, {1, 5, 0, 4, 3, 7, 2, 6}, {1, 5, 3, 7, 0, 4, 2, 6},
    {2, 0, 3, 1, 6, 4, 7, 5}, {2, 0, 6, 4, 3, 1, 7, 5}, {2, 3, 0, 1, 6, 7, 4, 5},
    {2, 3, 6, 7, 0, 1, 4, 5}, {2, 6, 0, 4, 3, 7, 1, 5}, {2, 6, 3, 7, 0, 4, 1, 5},
    {3, 1, 2, 0, 7, 5, 6, 4}, {3, 1, 7, 5, 2, 0, 6, 4}, {3, 2, 1, 0, 7, 6, 5, 4},
    {3, 2, 7, 6, 1, 0, 5, 4}, {3, 7, 1, 5, 2, 6, 0, 4}, {3, 7, 2, 6, 1, 5, 0, 4},
    {4, 0, 5, 1, 6, 2, 7, 3}, {4, 0, 6, 2, 5, 1, 7, 3}, {4, 5, 0, 1, 6, 7, 2, 3},
    {4, 5, 6, 7, 0, 1, 2, 3}, {4, 6, 0, 2, 5, 7, 1, 3}, {4, 6, 5, 7, 0, 2, 1, 3},
    {5, 1, 4, 0, 7, 3, 6, 2}, {5, 1, 7, 3, 4, 0, 6, 2}, {5, 4, 1, 0, 7, 6, 3, 2},
    {5, 4, 7, 6, 1, 0, 3, 2}, {5, 7, 1, 3, 4, 6, 0, 2}, {5, 7, 4, 6, 1, 3, 0, 2},
    {6, 2, 4, 0, 7, 3, 5, 1}, {6, 2, 7, 3, 4, 0, 5, 1}, {6, 4, 2, 0, 7, 5, 3, 1},
    {6, 4, 7, 5, 2, 0, 3, 1}, {6, 7, 2, 3, 4, 5, 0, 1}, {6, 7, 4, 5, 2, 3, 0, 1},
    {7, 3, 5, 1, 6, 2, 4, 0}, {7, 3, 6, 2, 5, 1, 4, 0}, {7, 5, 3, 1, 6, 4, 2, 0},
    {7, 5, 6, 4, 3, 1, 2, 0}, {7, 6, 3, 2, 5, 4, 1, 0}, {7, 6, 5, 4, 3, 2, 1, 0}
};

// Применить перестановку
std::array<Vector3d, 4> permute(const std::array<Vector3d, 4>& obj, int r) {
    return {
        obj[quad_symmetries[r][0]],
        obj[quad_symmetries[r][1]],
        obj[quad_symmetries[r][2]],
        obj[quad_symmetries[r][3]]
    };
}

// Применить перестановку
std::array<Vector3d, 8> permute(const std::array<Vector3d, 8>& obj, int r) {
    return {
        obj[cube_symmetries[r][0]],
        obj[cube_symmetries[r][1]],
        obj[cube_symmetries[r][2]],
        obj[cube_symmetries[r][3]],
        obj[cube_symmetries[r][4]],
        obj[cube_symmetries[r][5]],
        obj[cube_symmetries[r][6]],
        obj[cube_symmetries[r][7]]
    };
}

void subface_by_side_child_gen() {
    int gen_subface_by_side_child_2D[4][4] = {
        {Side2D::L[0], -1, Side2D::L[1], -1},
        {-1, Side2D::R[0], -1, Side2D::R[1]},
        {Side2D::B[0], Side2D::B[1], -1, -1},
        {-1, -1, Side2D::T[0], Side2D::T[1]}
    };
    // Проходим по subfaces
    for (Side2D s = 0; s < Side2D::n_subfaces(); ++s.val) {
        if (gen_subface_by_side_child_2D[s.simple()][s.child()] != s) {

            std::cout << s.simple() << " " << s.child() << " " << s << "\n";

            throw std::runtime_error("No way 2D");
        }
    }
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
        if (gen_subface_by_side_child_3D[s.simple()][s.child()] != s) {
            throw std::runtime_error("No way 3D");
        }
    }
    std::cout << "static constexpr int subface_by_side_child_3D[6][8] = {\n";
    for (Side3D side: Side3D::items()) {
        std::cout << "    {";
        for (int i = 0; i < 7; ++i) {
            std::cout << gen_subface_by_side_child_3D[side][i] << ", ";
        }
        std::cout << gen_subface_by_side_child_3D[side][7] << "},\n";
    }
    std::cout << "};\n";
}

void neib_child_gen_2D() {
    int gen_neib_child_by_symm_subface_2D[quad_symmetries.size()][Side2D::n_subfaces()];

    std::array<Vector3d, 4> UQ = {
        Vector3d{-0.5, -0.5, 0.0},
        Vector3d{+0.5, -0.5, 0.0},
        Vector3d{-0.5, +0.5, 0.0},
        Vector3d{+0.5, +0.5, 0.0}
    };

    SqQuad quad(UQ[0], UQ[1], UQ[2], UQ[3]);

    for (int r = 0; r < quad_symmetries.size(); ++r) {
        auto RQ = permute(UQ, r);

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

        std::array<std::array<Vector3d, CpC(2)>, Side2D::count()> child_c;
        for (Side2D side: Side2D::items()) {
            child_c[side] = {
                neibs[side](-0.5, -0.5),
                neibs[side](+0.5, -0.5),
                neibs[side](-0.5, +0.5),
                neibs[side](+0.5, +0.5)
            };
        }

        for (Side2D side = 0; side < Side2D::n_subfaces(); ++side.val) {
            Vector3d fc = 0.5 * (quad[side.cf()[0]] + quad[side.cf()[1]]);

            double min_dist = 100.0;
            int child_idx = -1;
            for (int i = 0; i < CpC(2); ++i) {
                double dist = (child_c[side.simple()][i] - fc).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    child_idx = i;
                }
            }
            gen_neib_child_by_symm_subface_2D[r][side] = child_idx;
        }
    }

    std::cout << "static constexpr int neib_child_by_symm_subface_2D[8][8] = {\n";
    for (int i = 0; i < quad_symmetries.size(); ++i) {
        std::cout << "    {";
        for (int j = 0; j < Side2D::n_subfaces() - 1; ++j) {
            std::cout << gen_neib_child_by_symm_subface_2D[i][j] << ", ";
        }
        std::cout << gen_neib_child_by_symm_subface_2D[i][Side2D::n_subfaces() - 1] << "},\n";
    }
    std::cout << "};\n";

}

void neib_child_gen_3D() {
    int gen_neib_child_by_symm_subface_3D[cube_symmetries.size()][Side3D::n_subfaces()];

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

    for (int r = 0; r < cube_symmetries.size(); ++r) {
        auto RQ = permute(UQ, r);

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

        std::array<std::array<Vector3d, CpC(3)>, Side3D::count()> child_c;
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
                quad[side.cf()[0]] + quad[side.cf()[1]] +
                quad[side.cf()[2]] + quad[side.cf()[3]]);

            double min_dist = 100.0;
            int child_idx = -1;
            for (int i = 0; i < CpC(3); ++i) {
                double dist = (child_c[side.simple()][i] - fc).norm();
                if (dist < min_dist) {
                    min_dist = dist;
                    child_idx = i;
                }
            }
            gen_neib_child_by_symm_subface_3D[r][side] = child_idx;
        }
    }

    std::cout << "static constexpr int neib_child_by_symm_subface_3D[48][24] = {\n";
    for (int i = 0; i < cube_symmetries.size(); ++i) {
        std::cout << "    {";
        for (int j = 0; j < Side3D::n_subfaces() - 1; ++j) {
            std::cout << gen_neib_child_by_symm_subface_3D[i][j] << ", ";
        }
        std::cout << gen_neib_child_by_symm_subface_3D[i][Side3D::n_subfaces() - 1] << "},\n";
    }
    std::cout << "};\n";
}

int main() {
    subface_by_side_child_gen();
    neib_child_gen_2D();
    neib_child_gen_3D();
    return 0;
}