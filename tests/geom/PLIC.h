// Header для тестов PLIC_2D и PLIC_3D
#pragma once

#include <iomanip>
#include <atomic>

#include <zephyr/math/funcs.h>

#include <zephyr/geom/geom.h>
#include <zephyr/geom/plic.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/io/pvd_file.h>

using namespace zephyr;
using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Cuboid;
using generator::Rectangle;

// val = max(val, count)
inline void update_max(std::atomic<size_t>& val, size_t count) {
    size_t current = val.load(std::memory_order_acquire);
    while (current < count &&
        !val.compare_exchange_weak(current, count,
            std::memory_order_release, std::memory_order_acquire)) { }
}

// Радиус и угол в декартовы координаты
inline Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

// Характеристическая функция (функция-индикатор)
using InFunction = std::function<bool(const Vector3d &)>;

// Пространственная функция
using SpFunction = std::function<double(const Vector3d &)>;
