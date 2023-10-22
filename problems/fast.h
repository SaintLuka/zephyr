#pragma once

#include <iostream>
#include <iomanip>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/geom/primitives/amr_face.h>
#include <zephyr/geom/vector.h>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using zephyr::mesh::generator::Strip;
using zephyr::mesh::generator::Rectangle;
using zephyr::mesh::generator::Cuboid;
using zephyr::mesh::generator::Sector;
using zephyr::mesh::Mesh;
using zephyr::mesh::Storage;
using zephyr::mesh::IFace;
using zephyr::mesh::ICell;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::geom::Matrix3d;

using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

using namespace zephyr::io;