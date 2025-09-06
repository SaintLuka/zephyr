
#include <zephyr/utils/mpi.h>

namespace zephyr::utils {

#ifdef ZEPHYR_MPI

static int g_size  = 1;
static int g_rank  = 0;
static int g_tasks = 1;

// Собирает уникальные имена процессоров/узлов
std::vector<std::string> proc_names();

void mpi_init(int& argc, char**& argv) {
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        return;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

    auto names = proc_names();
    std::string my_name = names[g_rank];

    g_tasks = 0;
    for (int r = 0; r < g_size; ++r) {
        if (names[r] == my_name) {
            g_tasks += 1;
        }
    }
}

mpi::handler::handler() {
    int argc = 0;
    char** argv = nullptr;

    mpi_init(argc, argv);
    default_types.init();
}

mpi::handler::handler(int &argc, char **&argv) {
    mpi_init(argc, argv);
    default_types.init();
}

mpi::handler::~handler() {
    default_types.free();

    int flag;
    MPI_Finalized(&flag);
    if (!flag) {
        MPI_Finalize();
    }
}

mpi::DefaultTypes mpi::default_types = {};

void mpi::DefaultTypes::init() {
    MPI_Type_contiguous(3, MPI_DOUBLE, &vec3);
    MPI_Type_commit(&vec3);

    MPI_Type_contiguous(4, MPI_INT, &int4);
    MPI_Type_commit(&int4);

    MPI_Type_contiguous(8, MPI_SHORT, &short8);
    MPI_Type_commit(&short8);

    byte4_types[0] = MPI_DATATYPE_NULL;
    for (int i = 1; i < byte4_count; ++i) {
        MPI_Type_contiguous(4 * i, MPI_BYTE, &byte4_types[i]);
        MPI_Type_commit(&byte4_types[i]);
    }
}

void mpi::DefaultTypes::free() {
    MPI_Type_free(&vec3);
    MPI_Type_free(&int4);
    MPI_Type_free(&short8);

    for (int i = 1; i < byte4_count; ++i) {
        MPI_Type_free(&byte4_types[i]);
    }
}

void mpi_free_type(MPI_Datatype& type) {
    MPI_Type_free(&type);
}

int mpi::rank() {
    return g_rank;
}

std::string mpi::srank() {
    return std::to_string(g_rank);
}

int mpi::size() {
    return g_size;
}

int mpi::n_tasks() {
    return g_tasks;
}

bool mpi::master() {
    return g_rank < 1;
}

bool mpi::single() {
    return g_size < 2;
}

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

// Собирает уникальные имена процессоров/узлов
std::vector<std::string> proc_names() {
    const int max_length = MPI_MAX_PROCESSOR_NAME;

    int length = 0;
    char name[max_length];

    MPI_Get_processor_name(name, &length);

    std::vector<int> lengths = mpi::all_gather(length);
    std::vector<int> offsets(lengths.size());
    offsets[0] = 0;
    for (int i = 1; i < g_size; ++i) {
        offsets[i] = offsets[i - 1] + lengths[i - 1];
    }

    int buff_size = offsets.back() + lengths.back();
    char *all_names = new char[buff_size + 1];

    MPI_Allgatherv(name, length, MPI_CHAR, all_names,
                   lengths.data(), offsets.data(), MPI_CHAR, MPI_COMM_WORLD);

    all_names[buff_size] = '\0';

    std::vector<std::string> names(g_size);
    for (int r = 0; r < g_size; ++r) {
        names[r] = std::string(all_names + offsets[r], lengths[r]);
    }

    //mpi::cout << "Process names:\n";
    //for (int r = 0; r < g_size; ++r) {
    //    mpi::cout << "  '" << names[r] << "'\n";
    //}

    return names;
}

mpi::master_stream mpi::cout = {};

namespace {

// Удаление MPI-типа, выполняет MPI_Type_free, а затем удаляет указатель на
// MPI_Datatype, добавляется к std::shared_ptr<MPI_Datatype>.
void free_deleter(MPI_Datatype *ptr) {
    MPI_Type_free(ptr);
    delete (ptr);
}
}

mpi::datatype::datatype(MPI_Datatype* dtype) {
    m_ptr = std::shared_ptr<MPI_Datatype>(dtype, free_deleter);
}

mpi::datatype::datatype(nullptr_t null) : m_ptr(nullptr) { }

mpi::datatype::datatype(MPI_Datatype dtype) : m_ptr(nullptr) {
    if (dtype != MPI_DATATYPE_NULL) {
        m_ptr = std::make_shared<MPI_Datatype>(dtype);
    }
}

mpi::datatype mpi::datatype::contiguous(int count, MPI_Datatype oldtype) {
    MPI_Datatype *newtype = new MPI_Datatype(MPI_DATATYPE_NULL);
    MPI_Type_contiguous(count, oldtype, newtype);
    MPI_Type_commit(newtype);
    return mpi::datatype(newtype);
}

mpi::datatype mpi::datatype::hvector(int count, int blocklength,
        MPI_Aint stride, MPI_Datatype oldtype) {
    MPI_Datatype *newtype = new MPI_Datatype(MPI_DATATYPE_NULL);
    MPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
    MPI_Type_commit(newtype);
    return mpi::datatype(newtype);
}

#else

std::ostream& mpi::cout = std::cout;

#endif

} // namespace zephyr::utils
