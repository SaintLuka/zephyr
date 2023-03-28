#include <control/mpi_wrapper.h>
#include <control/configuration.h>
#include <io/console.h>


int mpi::g_rank = 0;
int mpi::g_size = 1;

mpi::master_stream mpi::cout{};

bool exists_test(const string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// [1]
mpi::CellsTransaction::CellsTransaction(int target, int cells_count)
    : m_target(target)
    , m_cells_count(cells_count)
{ }

mpi::CellsTransaction::CellsTransaction(const CellsTransaction &transaction)
    : m_target(transaction.target())
    , m_cells_count(transaction.cells_count())
{ }

int mpi::CellsTransaction::target() const {
    return m_target;
}

int mpi::CellsTransaction::cells_count() const {
    return m_cells_count;
}
// [1]

/**
 * @brief Печать логотипа ASCII-ART
 */
void print_logo() {
    std::cout << R"( +======================================================================+ )" << std::endl;
    std::cout << R"( |      _   __  _   ___   _   __  _   _   __ __   ___   ___   _  _      | )" << std::endl;
    std::cout << R"( |     | | |  \| | | __| | | |  \| | | | |  V  | | __| /  _/ | || |     | )" << std::endl;
    std::cout << R"( |     | | | | ' | | _|  | | | | ' | | | | \_/ | | _|  `._`. | >< |     | )" << std::endl;
    std::cout << R"( |     |_| |_|\__| |_|   |_| |_|\__| |_| |_| |_| |___| |___/ |_||_|     | )" << std::endl;
    std::cout << R"( |                                                                      | )" << std::endl;
    std::cout << R"( +======================================================================+ )" << std::endl << std::endl;
}

void mpi::init(int argc, char **argv) {
#if IM_ENABLE_MPI
    if(argc == 0) {
        std::cerr << "Warning! You must create parallel controller before any object!" << std::endl;
        std::cerr << "Warning! Program will assume that no mpi running!" << std::endl;

        g_size = 1;
        g_rank = 0;
    }
    else {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &g_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
    }
#endif

    // Так как только тут определяется наш ранк - в этом же классе создаётся логгер
    auto filename = "process_" + to_string(rank()) + ".log";

    // Если такой уже есть - удалим
    if(exists_test(filename)) {
        std::ofstream ofs;
        ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
    }

    if (is_master()) {
        print_logo();
        std::cout << "  Parallel controller initialization\n";
        std::cout << "      Number of processes: " << g_size << "\n\n";
    }
}

void mpi::finilize() {
#if IM_ENABLE_MPI
    MPI_Finalize();
#endif
}

double mpi::min(double value) {
#if IM_ENABLE_MPI
    double glob_value = 0.0;
    MPI_Allreduce(&value, &glob_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

size_t mpi::min(size_t value) {
#if IM_ENABLE_MPI
    size_t glob_value = 0;
    MPI_Allreduce(&value, &glob_value, 1, MPI_SIZE_T, MPI_MIN, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

vector<double> mpi::min(const vector<double>& values) {
#if IM_ENABLE_MPI
    vector<double> glob_values(values.size());
    MPI_Allreduce(values.data(), glob_values.data(), static_cast<int>(values.size()), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return glob_values;
#else
    return values;
#endif
}

double mpi::max(double value) {
#if IM_ENABLE_MPI
    double glob_value = 0.0;
    MPI_Allreduce(&value, &glob_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

size_t mpi::max(size_t value) {
#if IM_ENABLE_MPI
    size_t glob_value = 0;
    MPI_Allreduce(&value, &glob_value, 1, MPI_SIZE_T, MPI_MAX, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

vector<double> mpi::max(const vector<double>& values) {
#if IM_ENABLE_MPI
    vector<double> glob_values(values.size());
    MPI_Allreduce(values.data(), glob_values.data(), static_cast<int>(values.size()), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return glob_values;
#else
    return values;
#endif
}

double mpi::sum(double value) {
#if IM_ENABLE_MPI
    double glob_value;
    MPI_Allreduce(&value, &glob_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

size_t mpi::sum(size_t value) {
#if IM_ENABLE_MPI
    size_t glob_value;
    MPI_Allreduce(&value, &glob_value, 1, MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    return glob_value;
#else
    return value;
#endif
}

vector<double> mpi::sum(const vector<double>& values) {
#if IM_ENABLE_MPI
    size_t size = values.size();
    size_t glob_size = 0;
    MPI_Allreduce(&size, &glob_size, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    if (glob_size != size) {
        throw runtime_error("Error in mpi::sum: Summation of different size vectors");
    }

    vector<double> glob_values(glob_size);
    MPI_Allreduce(values.data(), glob_values.data(), static_cast<int>(glob_size), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return glob_values;
#else
    return values;
#endif
}

mpi::loc_out<double> mpi::argmin(double value) {
#if IM_ENABLE_MPI
    mpi::loc_out<double> in(value, rank());
    mpi::loc_out<double> out(value, rank());
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    return out;
#else
    return mpi::loc_out<double>(value, rank());
#endif
}

mpi::loc_out<double> mpi::argmax(double value) {
#if IM_ENABLE_MPI
    mpi::loc_out<double> in(value, rank());
    mpi::loc_out<double> out(value, rank());
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    return out;
#else
    return mpi::loc_out<double>(value, rank());
#endif
}

vector<double> mpi::all_gather(double value) {
#if IM_ENABLE_MPI
    vector<double> values((size_t)g_size);
    MPI_Allgather(&value, 1, MPI_DOUBLE, values.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    return values;
#else
    return std::VECTOR<double>(1, value);
#endif
}

vector<size_t> mpi::all_gather(size_t value) {
#if IM_ENABLE_MPI
    vector<size_t> values((size_t)g_size);
    MPI_Allgather(&value, 1, MPI_SIZE_T, values.data(), 1, MPI_SIZE_T, MPI_COMM_WORLD);
    return values;
#else
    return std::VECTOR<size_t>(1, value);
#endif
}

vector<double> mpi::all_gather(const vector<double> &values) {
#if IM_ENABLE_MPI
    size_t all_sizes;
    size_t size = values.size();
    MPI_Allreduce(&size, &all_sizes, 1, MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

    vector<double> all_values(all_sizes);
    MPI_Allgather(values.data(), (int)size, MPI_DOUBLE, all_values.data(), (int)size, MPI_DOUBLE, MPI_COMM_WORLD);
    return all_values;
#else
    return values;
#endif
}

vector<size_t> mpi::all_gather(const array<size_t, 3> &values) {
#if IM_ENABLE_MPI
    size_t all_sizes;
    size_t size = values.size();
    MPI_Allreduce(&size, &all_sizes, 1, MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

    vector<size_t> all_values(all_sizes);
    MPI_Allgather(values.data(), (int)size, MPI_SIZE_T, all_values.data(), (int)size, MPI_SIZE_T, MPI_COMM_WORLD);
    return all_values;
#else
    vector<size_t> valuess = {values[0], values[1], values[2]};
    return valuess;
#endif
}

int mpi::rank() {
    return g_rank;
}

unsigned int mpi::urank() {
    return static_cast<unsigned int>(g_rank);
}

const string& mpi::srank() {
    static string s_rank = to_string(mpi::g_rank);
    return s_rank;
}

int mpi::size() {
    return g_size;
}

unsigned int mpi::usize() {
    return static_cast<unsigned int>(g_size);
}

bool mpi::is_master() {
    return g_rank == 0;
}

void mpi::barrier() {
#if IM_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

