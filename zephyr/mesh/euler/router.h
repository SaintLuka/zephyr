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
    FACE_BEG,
    NODE_BEG,

    // Данные граней
    ADJ_RANK,
    ADJ_INDEX,
    ADJ_ALIEN,
    ADJ_BASIC,
    BOUNDARY,
    NORMAL,
    FACE_CENTER,
    AREA,
    AREA_ALT,
    FACE_VERTS,

    // Данные вершин
    VERTICES
};

inline std::string to_string(MpiTag tag) {
    return std::to_string(static_cast<int>(tag));
}

#ifdef ZEPHYR_MPI

/// @brief Массив MPI-запросов isend/irecv
class Requests {
public:
    /// @brief Массив "нулевых" запросов
    explicit Requests()
        : m_requests(utils::mpi::size(), MPI_REQUEST_NULL) { }

    /// @brief Массив "нулевых" запросов
    explicit Requests(int size)
        : m_requests(size, MPI_REQUEST_NULL) { }

    /// @brief Получить запрос по номеру ранга
    MPI_Request& operator[](int r) { return m_requests[r]; }

    /// @brief Дождаться завершения всех запросов
    void wait() const;

private:
    std::vector<MPI_Request> m_requests;
};

/// @brief Управляет обменными операциями
class Router {
public:
    /// @brief По умолчанию, size = mpi::size()
    /// Массивы инициализируются нулями
    Router();

    /// @brief Конструктор без инициализации массивов
    //explicit Router(int size);

    /// @brief Установить число элементов для отправки
    void set_send_count(const std::vector<index_t>& send_count);

    /// @brief Установить число элементов для получения
    void set_recv_count(const std::vector<index_t>& recv_count);

    /// @brief Собрать пересылки
    void fill_partial();

    /// @brief Собрать полную матрицу пересылок
    void fill_complete();

    /// @brief Число процессов (== mpi::size())
    int size() const { return m_size; }

    /// @brief Есть полная матрица пересылок?
    bool complete() const { return !m_send_recv.empty(); }


    /// @brief Количество пересылок с i-го процесса на j-ый
    index_t get(int i, int j) const;

    /// @brief Количество пересылок с i-го процесса на j-ый
    index_t operator()(int i, int j) const;

    /// @brief Необходимый размер буфера для отправки сообещений
    index_t send_buffer_size() const;

    /// @brief Необходимый размер буфера для получения сообещений
    index_t recv_buffer_size() const;

    const std::vector<index_t>& send_count() const { return m_send_count; }
    const std::vector<index_t>& recv_count() const { return m_recv_count; }

    const std::vector<index_t>& send_offset() const { return m_send_offset; }
    const std::vector<index_t>& recv_offset() const { return m_recv_offset; }

    index_t send_count(int r) const { return m_send_count[r]; }
    index_t recv_count(int r) const { return m_recv_count[r]; }

    index_t send_offset(int r) const { return m_send_offset[r]; }
    index_t recv_offset(int r) const { return m_recv_offset[r]; }

    /// @brief Индексы из массива m_border_indices
    range<index_t> send_indices(int r) const {
        return {m_send_offset[r], m_send_offset[r] + m_send_count[r]};
    }
    /// @brief Индексы из массива aliens при получении
    range<index_t> recv_indices(int r) const {
        return {m_recv_offset[r], m_recv_offset[r] + m_recv_count[r]};
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

    Requests isend(const utils::Buffer& src, MpiTag tag);


    /// @brief Асинхронное получение
    template<typename T>
    Requests irecv(T* dst, MpiTag tag, MPI_Datatype dtype);

    template<typename T>
    Requests irecv(T* dst, MpiTag tag);

    template<typename T>
    Requests irecv(std::vector<T>& dst, MpiTag tag);

    template<typename T>
    Requests irecv(std::vector<T>& dst, MpiTag tag, MPI_Datatype dtype);

    Requests irecv(utils::Buffer& dst, MpiTag tag);

protected:

    /// @brief Вывести полную матрицу пересылок в консоль
    void print_partial() const;

    /// @brief Вывести полную матрицу пересылок в консоль
    void print_complete() const;


    int m_size;  ///< Число процессов (== mpi::size())

    std::vector<index_t> m_send_count;   ///< Число элементов на отправку
    std::vector<index_t> m_recv_count;   ///< Число элементов на получение
    std::vector<index_t> m_send_offset;  ///< Смещения в массиве на отправку
    std::vector<index_t> m_recv_offset;  ///< Смещения в массиве на получение

    /// @brief Полная матрица пересылок (опционально)
    std::vector<index_t> m_send_recv;
};


template<typename T>
Requests Router::isend(const T* src, MpiTag tag, MPI_Datatype dtype) {
    Requests send_req(m_size);
    for (int r = 0; r < m_size; ++r) {
        if (m_send_count[r] > 0) {
            MPI_Isend(src + m_send_offset[r], m_send_count[r],
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
    Requests recv_req(m_size);
    for (int r = 0; r < m_size; ++r) {
        if (m_recv_count[r] > 0) {
            MPI_Irecv(dst + m_recv_offset[r], m_recv_count[r],
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