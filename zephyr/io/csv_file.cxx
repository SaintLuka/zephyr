#include <fstream>
#include <iomanip>

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/io/csv_file.h>


namespace zephyr { namespace io {

using namespace zephyr::geom;

/// ===========================================================================
///             Реализация записи в виде набора статических функций
/// ===========================================================================

inline size_t count_cells(AmrStorage &cells, const Filter &filter) {
    if (filter.is_trivial()) {
        return cells.size();
    }
    size_t count = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            ++count;
        }
    }
    return count;
}

void write_name(const Variable& var, std::ofstream& file) {
    if (var.is_scalar()) {
        file << var.name();
    }
    else {
        for (int k = 0; k < var.n_components() - 1; ++k) {
            file << var.name() << "[" << k << "], ";
        }
        file << var.name() << "[" << var.n_components() - 1 << "]";
    }
}

void write_variable(const Variable& var, AmrStorage::Item& cell, std::ofstream& file) {
    if (!var.is_scalar()) {
        throw std::runtime_error("CsvFile doesn't support vector variables");
    }
    double out;

    var.write(cell, &out);

    switch (var.type()) {
        case VtkType::UInt8:   file << ((uint8_t*) (&out))[0]; break;
        case VtkType::UInt16:  file << ((uint16_t*)(&out))[0]; break;
        case VtkType::UInt32:  file << ((uint32_t*)(&out))[0]; break;
        case VtkType::UInt64:  file << ((uint64_t*)(&out))[0]; break;
        case VtkType::Int8:    file << ((int8_t*)  (&out))[0]; break;
        case VtkType::Int16:   file << ((int16_t*) (&out))[0]; break;
        case VtkType::Int32:   file << ((int32_t*) (&out))[0]; break;
        case VtkType::Int64:   file << ((int64_t*) (&out))[0]; break;
        case VtkType::Float32: file << ((float*)   (&out))[0]; break;
        case VtkType::Float64: file << ((double*)  (&out))[0]; break;
        default:
            throw std::runtime_error("CsvFile::write_variable: has no type");
    }
}

CsvFile::CsvFile(
        const std::string &filename, int precision,
        const Variables &variables) :
    filename(filename), precision(precision),
    variables(variables), filter() {
}

void CsvFile::save(AmrStorage &cells) {
    save(filename, cells, precision, variables, filter);
}

void CsvFile::save(
    const std::string &filename, AmrStorage &cells, int precision,
    const Variables &variables, const Filter &filter
) {
    size_t n_cells = count_cells(cells, filter);
    if (n_cells < 1) {
        return;
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    int n_variables = int(variables.size()) - 1;

    file << std::scientific << std::setprecision(precision);

    file << "# x, y, ";
    for(int i = 0; i < n_variables - 1; ++i) {
        write_name(variables[i], file);
        file << ", ";
    }
    if (n_variables > 0) {
        write_name(variables[n_variables - 1], file);
        file << "\n";
    }

    for (auto& cell: cells) {
        if (filter(cell)) {
            file << cell.center.x() << ", " << cell.center.y() << ", ";
            for(int i = 0; i < n_variables - 1; ++i) {
                write_variable(variables[i], cell, file);
                file << ", ";
            }
            if (n_variables > 0) {
                write_variable(variables[n_variables - 1], cell, file);
            }
            file << "\n";
        }
    }

    file.close();
}

} // io
} // zephyr