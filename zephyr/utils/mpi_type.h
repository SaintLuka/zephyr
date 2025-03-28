#pragma once

#include <zephyr/configuration.h>

#include <mpi.h>
#include <tuple>
#include <type_traits>

namespace zephyr::utils {

template<class T>
static MPI_Datatype mpi_type();

template <> MPI_Datatype mpi_type<int>() { return MPI_INT; }
template <> MPI_Datatype mpi_type<unsigned char>() { return MPI_BYTE; }
template <> MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
template <> MPI_Datatype mpi_type<float>() { return MPI_FLOAT; }
template <> MPI_Datatype mpi_type<long>() { return MPI_LONG; }
template <> MPI_Datatype mpi_type<short>() { return MPI_SHORT; }
        
class MPITypeManager {
public:
    static std::vector<MPI_Datatype>& get_types(){
        static std::vector<MPI_Datatype> types;
        return types;
    }

    static void commit_types(){
        for(MPI_Datatype& type : get_types()){
            MPI_Type_commit(&type);
        }
    }

    static void free_types(){
        for(MPI_Datatype& type : get_types()){
            MPI_Type_free(&type);
        }
    }

    static void register_type(MPI_Datatype type) {
        get_types().push_back(type);
    }
};

template<typename T, typename... Fields, size_t... Indices> 
MPI_Datatype _mpi_register_type(std::index_sequence<Indices...>) {
    constexpr int count = sizeof...(Fields);
    int block_lengths[count] = {1};
    MPI_Aint displacements[count];
    MPI_Datatype types[count] = {mpi_type<Fields>()...};

    std::tuple<Fields...> tmp;
    MPI_Aint base_address;
    MPI_Get_address(&tmp, &base_address);
    ((MPI_Get_address(&std::get<Indices>(tmp), &displacements[Indices]), displacements[Indices] -= base_address), ...);

    MPI_Datatype new_type;
    MPI_Type_create_struct(count, block_lengths, displacements, types, &new_type);
    MPI_Type_commit(&new_type);
    MPITypeManager::register_type(new_type);

    return new_type;
}

template<typename T, typename... Fields>
MPI_Datatype mpi_register_type() {
    static MPI_Datatype type = _mpi_register_type<T, Fields...>(std::index_sequence_for<Fields...>{});
    return type;
}

}