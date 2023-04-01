#pragma once

#include <vector>

#include <zephyr/configuration.h>

#ifdef ZEPHYR_ENABLE_MPI
#include <mpi.h>
#endif


namespace zephyr { namespace utils {

/// @brief Надо придумать, как это изящнее реализовать
template <class T>
MPI_Datatype mpi_type() {
    return MPI_INT;
}

/// @brief Статический класс. Упрощенная обертка для mpi.
/// @details При выключенном mpi все функции продолжают работать,
/// так что их можно не удалять из однопроцессорного кода.
class mpi {
public:

    /// @defgroup Базовые функции
    /// @{

    /// @brief Инициализация mpi
    static void init();

    /// @brief Завершение
    static void finalize();

    /// @brief Ранг процесса
    static int rank();

    /// @brief Ранг процесса в виде строки
    static std::string srank();

    /// @brief Число процессов
    static int size();

    /// @brief Это мастер процесс?
    static bool is_master();

    /// @brief Приложение однопроцессорное?
    static bool is_single();

#ifdef ZEPHYR_ENABLE_MPI
    /// @brief Сеть коммуникаций
    inline static MPI_Comm comm() { return MPI_COMM_WORLD; };
#endif

    /// @brief Блокирует все процессы данной сети, пока они не достигнут
    /// данного вызова.
    static void barrier();

    /// @brief Выполнить функцию последоватьно (!) на каждом процессе
    /// @details Помогает при отладке
    /// @param F Целевая функция
    template <class F>
    static void for_each(F&&);

    /// @brief Прикольный интерфейс, чтобы можно было писать
    /// mpi::cout << ... и сообщения выводил только мастер-процесс
    struct master_stream {
        template <class T>
        std::ostream& operator<<(const T& val) {
            if (mpi::is_master()) {
                std::cout << val;
            }
            return std::cout;
        }
    };

    /// @brief Просто пиши mpi::cout <<
    static master_stream cout;

    /// @}

    /// @defgroup Коллективные reduce операции
    /// @{

    /// @brief Коллективная операция. Минимальное значение по всем процессам.
    template <class T>
    static T min(const T& value);

    /// @brief Коллективная операция. Покомпонентный минимум для каждого
    /// элемента вектора серди всех процессов сети.
    template <class T>
    static std::vector<T> min(const std::vector<T>& values);

    /// @brief Коллективная операция. Максимальное значение по всем процессам.
    template <class T>
    static T max(const T& value);

    /// @brief Коллективная операция. Покомпонентный максимум для каждого
    /// элемента вектора серди всех процессов сети.
    template <class T>
    static T max(const std::vector<T>& values);

    /// @brief Коллективная операция. Сумма величин со всех процессов сети.
    template <class T>
    static T sum(const T& value);

    /// @brief Коллективная операция. Покомпонентная сумма элементов
    // векторов со всех процессов сети.
    template <class T>
    static std::vector<T> sum(const std::vector<T>& values);

    /// @}

    /// @defgroup Коллективные обменные операции
    /// @{

    /// @brief Отправляет сообщение с "корневого" процесса всем процессам.
    /// @param root Ранг "корневого" процесса.
    /// @param value Величина шаблонного типа по ссылке, переменная также
    /// принимает пересланное значение.
    template <class T>
    static void broadcast(int root, T& value);

    /// @brief Коллективная операция. Собирает величину value со всех
    /// процессов, записывает в массив и раздает этот массив всем процессам.
    template <class T>
    static std::vector<T> all_gather(const T& value);

    /// @brief Коллективная операция. Собирает величину value со всех процессов,
    /// записывает в массив и раздает этот массив всем процессам сети.
    template <class T>
    static void all_gather(const T& value, std::vector<T>& values);

    /// @brief Коллективная операция. Размеры буфферов send и recv совпадают
    /// и равны числу процессов. Каждая величина из send отправляется своему
    /// процессу, каждая величина в recv получается с определенного процесса.
    template <class T>
    static void all_to_all(const std::vector<T>& send, std::vector<T>& recv);

    /// @}
};

template <class F>
void mpi::for_each(F&& f) {
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == mpi::rank()) {
            f();
            std::cout << std::flush;
        }
        mpi::barrier();
    }
}

template <class T>
T mpi::min(const T& value) {
#ifdef ZEPHYR_ENABLE_MPI
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MIN, comm());
    return g_value;
#else
    return value;
#endif
}

template <class T>
std::vector<T> mpi::min(const std::vector<T>& values) {
#ifdef ZEPHYR_ENABLE_MPI
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MIN, comm());
    return g_values;
#else
    return values;
#endif
}

template <class T>
T mpi::max(const T& value) {
#ifdef ZEPHYR_ENABLE_MPI
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MAX, comm());
    return g_value;
#else
    return value;
#endif
}

template <class T>
T mpi::max(const std::vector<T>& values) {
#ifdef ZEPHYR_ENABLE_MPI
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MAX, comm());
    return g_values;
#else
    return values;
#endif
}

template <class T>
T mpi::sum(const T& value) {
#ifdef ZEPHYR_ENABLE_MPI
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_SUM, comm());
    return g_value;
#else
    return value;
#endif
}

template <class T>
std::vector<T> mpi::sum(const std::vector<T>& values) {
#ifdef ZEPHYR_ENABLE_MPI
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_SUM, comm());
    return g_values;
#else
    return values;
#endif
}

template <class T>
void mpi::broadcast(int root, T& value) {
#ifdef ZEPHYR_ENABLE_MPI
    MPI_Bcast((void*)&value, sizeof(value), MPI_CHAR, root, comm());
#endif
}

template <class T>
std::vector<T> mpi::all_gather(const T& value) {
#ifdef ZEPHYR_ENABLE_MPI
    std::vector<T> values(size());
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), comm());
    return values;
#else
    return { value };
#endif
}

template <class T>
void mpi::all_gather(const T& value, std::vector<T>& values) {
#ifdef ZEPHYR_ENABLE_MPI
    values.resize(size());
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), comm());
#else
    values = { value };
#endif
}

template <class T>
void mpi::all_to_all(const std::vector<T>& send, std::vector<T>& recv) {
#ifdef ZEPHYR_ENABLE_MPI
    recv.resize(size());
    MPI_Alltoall(send.data(), 1, mpi_type<T>(), recv.data(), 1, mpi_type<T>(), comm());
#else
    recv = send;
#endif
}

} // utils
} // zephyr