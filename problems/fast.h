#pragma once

#include <iostream>
#include <iomanip>

#include <zephyr/geom/face.h>
#include <zephyr/geom/vector.h>

#include <zephyr/mesh/generator/rectangle.h>
#include <zephyr/mesh/generator/sector.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/io/pvd_file.h>

using zephyr::mesh::generator::Rectangle;
using zephyr::mesh::generator::Sector;
using zephyr::mesh::Mesh;
using zephyr::mesh::Storage;
using zephyr::geom::FaceFlag;
using zephyr::geom::Vector3d;
using zephyr::geom::Matrix3d;

using namespace zephyr::io;