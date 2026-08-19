#pragma once

#include <vector>
#include <future>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/buffer.h>
#include <zephyr/utils/range.h>
#include <zephyr/mesh/index.h>

namespace zephyr::mesh {

/// @brief MPI-тэги для корректных пересылок
enum class MpiTag : int {
    NONE = 0,

    // Данные ячеек
    RANK,
    NEXT,
    INDEX,
    FLAG,
    LEVEL,
    B_IDX,
    Z_IDX,
    CENTER,
    VOLUME,
    VOLUME_ALT,

    // Данные граней
    FACE_BEG,
    ADJ_RANK,
    ADJ_INDEX,
    ADJ_GHOST,
    ADJ_BASIC,
    ADJ_ROTATION,
    BOUNDARY,
    NORMAL,
    FACE_CENTER,
    AREA,
    AREA_ALT,
    FACE_VERTS,

    // Данные вершин
    VERT_BEG,
    VERT_COORD,
    VERT_INDEX,
    VERT_GHOST
};

inline std::string to_string(MpiTag tag) {
    return std::to_string(static_cast<int>(tag));
}

#ifdef ZEPHYR_MPI

/// @brief Массив MPI-запросов isend/irecv
class Requests {
public:
    /// @brief Массив "нулевых" запросов
    explicit Requests(int size);

    /// @brief Получить запрос по номеру ранга
    MPI_Request& operator[](int r) { return requests_[r]; }

    /// @brief Дождаться завершения всех запросов
    void wait() const;

private:
    /// @brief Число MPI-процессов
    int size_;

    /// @brief Массив запросов для каждого процесса
    /// (фактически запрет копирования)
    std::unique_ptr<MPI_Request[]> requests_;
};

class RequestsList {
public:
    RequestsList() = default;

    void reserve(int size);

    void operator+=(Requests&& requests);

    void wait() const;

private:
    std::vector<Requests> requests_;
};

/// @brief Управляет обменными операциями
class Router {
public:
    /// @brief По умолчанию, size = mpi::size()
    /// Массивы инициализируются нулями
    Router();

    /// @brief Установить число элементов для отправки
    void set_send_count(const std::vector<index_t>& send_count);

    /// @brief Установить число элементов для получения
    void set_recv_count(const std::vector<index_t>& recv_count);

    /// @brief Собрать пересылки
    void fill_partial();

    /// @brief Собрать полную матрицу пересылок
    void fill_complete();

    /// @brief Число процессов (== mpi::size())
    int size() const { return size_; }

    /// @brief Есть полная матрица пересылок?
    bool complete() const { return !send_recv_.empty(); }


    /// @brief Количество пересылок с i-го процесса на j-ый
    index_t get(int i, int j) const;

    /// @brief Количество пересылок с i-го процесса на j-ый
    index_t operator()(int i, int j) const;

    /// @brief Необходимый размер буфера для отправки сообщений
    index_t send_buffer_size() const;

    /// @brief Необходимый размер буфера для получения сообщений
    index_t recv_buffer_size() const;

    const std::vector<index_t>& send_count() const { return send_count_; }
    const std::vector<index_t>& recv_count() const { return recv_count_; }

    const std::vector<index_t>& send_offset() const { return send_offset_; }
    const std::vector<index_t>& recv_offset() const { return recv_offset_; }

    index_t send_count(int r) const { return send_count_[r]; }
    index_t recv_count(int r) const { return recv_count_[r]; }

    index_t send_offset(int r) const { return send_offset_[r]; }
    index_t recv_offset(int r) const { return recv_offset_[r]; }

    /// @brief Индексы из массива border_indices_
    range_t<index_t> send_indices(int r) const {
        return range(send_offset_[r], send_offset_[r] + send_count_[r]);
    }

    /// @brief Индексы из массива ghosts при получении
    range_t<index_t> recv_indices(int r) const {
        return range(recv_offset_[r], recv_offset_[r] + recv_count_[r]);
    }

    /// @brief Вывести информацию о пересылках
    void print() const;


    /// @brief Асинхронная отправка
    template<typename T>
    Requests isend(const T* src, MpiTag tag, MPI_Datatype dtype);

    template<typename T>
    Requests isend(const T* src, MpiTag tag);

    template<typename T>
    Requests isend(const std::vector<T>& src, MpiTag tag);

    template<typename T>
    Requests isend(const std::vector<T>& src, MpiTag tag, MPI_Datatype dtype);

    Requests isend(const utils::Buffer& src, MpiTag tag) const;


    /// @brief Асинхронное получение
    template<typename T>
    Requests irecv(T* dst, MpiTag tag, MPI_Datatype dtype);

    template<typename T>
    Requests irecv(T* dst, MpiTag tag);

    template<typename T>
    Requests irecv(std::vector<T>& dst, MpiTag tag);

    template<typename T>
    Requests irecv(std::vector<T>& dst, MpiTag tag, MPI_Datatype dtype);

    Requests irecv(utils::Buffer& dst, MpiTag tag) const;

protected:

    /// @brief Вывести полную матрицу пересылок в консоль
    void print_partial() const;

    /// @brief Вывести полную матрицу пересылок в консоль
    void print_complete() const;


    int size_;  ///< Число процессов (== mpi::size())

    std::vector<index_t> send_count_;   ///< Число элементов на отправку
    std::vector<index_t> recv_count_;   ///< Число элементов на получение
    std::vector<index_t> send_offset_;  ///< Смещения в массиве на отправку
    std::vector<index_t> recv_offset_;  ///< Смещения в массиве на получение

    /// @brief Полная матрица пересылок (опционально)
    std::vector<index_t> send_recv_;
};


template<typename T>
Requests Router::isend(const T* src, MpiTag tag, MPI_Datatype dtype) {
    Requests send_req(size_);
    for (int r = 0; r < size_; ++r) {
        if (send_count_[r] > 0) {
            MPI_Isend(src + send_offset_[r], send_count_[r],
                      dtype, r, int(tag), utils::mpi::comm(), &send_req[r]);
        }
    }
    return send_req;
}

template<typename T>
Requests Router::isend(const T* src, MpiTag tag) {
    return isend(src, tag, utils::mpi::type<T>());
}

template<typename T>
Requests Router::isend(const std::vector<T>& src, MpiTag tag) {
    return isend(src, tag, utils::mpi::type<T>());
}

template<typename T>
Requests Router::isend(const std::vector<T>& src, MpiTag tag, MPI_Datatype dtype) {
    return isend(src.data(), tag, dtype);
}

template<typename T>
Requests Router::irecv(T* dst, MpiTag tag, MPI_Datatype dtype) {
    Requests recv_req(size_);
    for (int r = 0; r < size_; ++r) {
        if (recv_count_[r] > 0) {
            MPI_Irecv(dst + recv_offset_[r], recv_count_[r],
                      dtype, r, int(tag), utils::mpi::comm(), &recv_req[r]);
        }
    }
    return recv_req;
}

template<typename T>
Requests Router::irecv(T* dst, MpiTag tag) {
    return irecv(dst, tag, utils::mpi::type<T>());
}

template<typename T>
Requests Router::irecv(std::vector<T>& dst, MpiTag tag) {
    return irecv(dst, tag, utils::mpi::type<T>());
}

template<typename T>
Requests Router::irecv(std::vector<T>& dst, MpiTag tag, MPI_Datatype dtype) {
    return irecv(dst.data(), tag, dtype);
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh