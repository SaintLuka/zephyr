#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <zephyr/configuration.h>

#ifdef ZEPHYR_MPI

#include <mpi.h>

// Нормальная параллельная версия обертки mpi

namespace zephyr::utils {

/// @brief Статический класс. Упрощенная обертка для mpi.
class mpi {
public:

    /// @{
    /// @name Некоторое пояснение к обычной группе? Базовые функции, да как это работает???

    /// @brief Инициализация mpi
    static void init();

    /// @brief Завершение
    static void finalize();

    // @brief Rfrjuj [thf Ранг процесса
    static int rank();

    /// @brief Ранг процесса в виде строки
    static std::string srank();

    /// @brief Число процессов
    static int size();
    
    /// @brief Число задач на узел
    static int n_tasks();

    /// @brief Это мастер процесс?
    static bool master();

    /// @brief Приложение однопроцессорное?
    static bool single();

    /// @brief Сеть коммуникаций
    inline static constexpr MPI_Comm comm() { return MPI_COMM_WORLD; };

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
        master_stream& operator<<(const T& val) {
            if (mpi::master()) {
                std::cout << val;
            }
            return *this;
        }
    };

    /// @brief Просто пиши mpi::cout <<
    static master_stream cout;

    /// @}

    /// @{
    /// @name Коллективные reduce операции

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
    static std::vector<T> max(const std::vector<T>& values);

    /// @brief Коллективная операция. Сумма величин со всех процессов сети.
    template <class T>
    static T sum(const T& value);

    /// @brief Коллективная операция. Покомпонентная сумма элементов
    // векторов со всех процессов сети.
    template <class T>
    static std::vector<T> sum(const std::vector<T>& values);

    /// @brief Номер процесса и значение
    template <class T>
    struct value_rank { T value; int rank; };

    /// @brief Минимальное значение и номер процесса, на котором достигается
    template <class T>
    static value_rank<T> argmin(const T& value);

    /// @brief Максимальное значение и номер процесса, на котором достигается
    template <class T>
    static value_rank<T> argmax(const T& value);

    /// @}

    /// @{
    /// @name Коллективные обменные операции

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

    /// @brief Коллективная операция. Собирает вектор величин value со всех процессов,
    /// записывает в один общий массив и раздает этот массив всем процессам сети.
    template <class T>
    static std::vector<T> all_gather_vectors(const std::vector<T>& value);

    /// @brief Коллективная операция. Размеры буфферов send и recv совпадают
    /// и равны числу процессов. Каждая величина из send отправляется своему
    /// процессу, каждая величина в recv получается с определенного процесса.
    template <class T>
    static void all_to_all(const std::vector<T>& send, std::vector<T>& recv);

    /// @}
};

template<class T>
static MPI_Datatype mpi_type();

template <> MPI_Datatype mpi_type<int>() { return MPI_INT; }
template <> MPI_Datatype mpi_type<unsigned char>() { return MPI_BYTE; }
template <> MPI_Datatype mpi_type<char>() { return MPI_BYTE; }
template <> MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
template <> MPI_Datatype mpi_type<float>() { return MPI_FLOAT; }
template <> MPI_Datatype mpi_type<long>() { return MPI_LONG; }
template <> MPI_Datatype mpi_type<short>() { return MPI_SHORT; }

template<class T>
MPI_Datatype mpi_register_contiguous_type(int count) {
    MPI_Datatype out;
    MPI_Type_contiguous(
        count,
        mpi_type<T>(),
        &out
    );
    MPI_Type_commit(&out);
    return out;
}

void mpi_free_type(MPI_Datatype& type);

template <class F>
void mpi::for_each(F&& f) {
    if (mpi::single()) {
        f(); std::cout << std::flush;
        return;
    }
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
    if (mpi::single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MIN, comm());
    return g_value;
}

template <class T>
std::vector<T> mpi::min(const std::vector<T>& values) {
    if (mpi::single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MIN, comm());
    return g_values;
}

template <class T>
T mpi::max(const T& value) {
    if (mpi::single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MAX, comm());
    return g_value;
}

template <class T>
std::vector<T> mpi::max(const std::vector<T>& values) {
    if (mpi::single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MAX, comm());
    return g_values;
}

template <class T>
T mpi::sum(const T& value) {
    if (mpi::single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_SUM, comm());
    return g_value;
}

template <class T>
std::vector<T> mpi::sum(const std::vector<T>& values) {
    if (mpi::single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_SUM, comm());
    return g_values;
}

template <class T>
mpi::value_rank<T> mpi::argmin(const T& value) {
    static_assert(std::is_same<T, double>::value, "Only doubles");
    mpi::value_rank<T> in {.value=value, .rank=rank()};
    mpi::value_rank<T> out{.value=value, .rank=rank()};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    return out;
}

template <class T>
mpi::value_rank<T> mpi::argmax(const T& value) {
    static_assert(std::is_same<T, double>::value, "Only doubles");
    mpi::value_rank<T> in {.value=value, .rank=rank()};
    mpi::value_rank<T> out{.value=value, .rank=rank()};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    return out;
}

template <class T>
void mpi::broadcast(int root, T& value) {
    MPI_Bcast((void*)&value, sizeof(value), MPI_CHAR, root, comm());
}

template <class T>
std::vector<T> mpi::all_gather(const T& value) {
    std::vector<T> values(size());
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), comm());
    return values;
}

template <class T>
void mpi::all_gather(const T& value, std::vector<T>& values) {
    values.resize(size());
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), comm());
}

template <class T>
std::vector<T> mpi::all_gather_vectors(const std::vector<T>& value) {
    std::vector<T> values(size() * value.size());
    MPI_Allgather(value.data(), value.size(), mpi_type<T>(), values.data(), value.size(), mpi_type<T>(), comm());
    return values;
}

template <class T>
void mpi::all_to_all(const std::vector<T>& send, std::vector<T>& recv) {
    recv.resize(size());
    MPI_Alltoall(send.data(), 1, mpi_type<T>(), recv.data(), 1, mpi_type<T>(), comm());
}

} // namespace zephyr::utils

#else

// Однопроцессорная заглушка при компиляции без mpi

namespace zephyr::utils {

/// @brief Статический класс. Заглушка для mpi.
class mpi {
public:
    /// @{
    /// @name Базовые функции MPI

    /// @brief Инициализация mpi
    static void init() { };

    /// @brief Завершение
    static void finalize() { };

    /// @brief Ранг процесса
    static constexpr int rank() { return 0; }

    /// @brief Ранг процесса в виде строки
    static std::string srank() { return "0"; };

    /// @brief Число процессов
    static constexpr int size() { return 1; };
    
    /// @brief Число задач на узел
    static constexpr int n_tasks() { return 1; }

    /// @brief Это мастер процесс?
    static constexpr bool master() { return true; }

    /// @brief Приложение однопроцессорное?
    static constexpr bool single() { return true; }

    /// @brief Сеть коммуникаций
    inline static constexpr bool comm() { return false; };

    /// @brief Блокирует все процессы данной сети, пока они не достигнут
    /// данного вызова.
    static void barrier() { };

    /// @brief Выполнить функцию последоватьно (!) на каждом процессе
    /// @details Помогает при отладке
    /// @param F Целевая функция
    template <class F>
    static void for_each(F&& f) { f(); }

    /// @brief Просто пиши mpi::cout <<
    static std::ostream& cout;

    /// @}

    /// @{
    /// @name Коллективные reduce операции

    /// @brief Коллективная операция. Минимальное значение по всем процессам.
    template <class T>
    static T min(const T& value) { return value; }

    /// @brief Коллективная операция. Покомпонентный минимум для каждого
    /// элемента вектора серди всех процессов сети.
    template <class T>
    static std::vector<T> min(const std::vector<T>& values) {
        return values;
    }

    /// @brief Коллективная операция. Максимальное значение по всем процессам.
    template <class T>
    static T max(const T& value) { return value; }

    /// @brief Коллективная операция. Покомпонентный максимум для каждого
    /// элемента вектора серди всех процессов сети.
    template <class T>
    static T max(const std::vector<T>& values) { return values; }

    /// @brief Коллективная операция. Сумма величин со всех процессов сети.
    template <class T>
    static T sum(const T& value) { return value; }

    /// @brief Коллективная операция. Покомпонентная сумма элементов
    // векторов со всех процессов сети.
    template <class T>
    static std::vector<T> sum(const std::vector<T>& values) { return values; }

    /// @brief Номер процесса и значение
    template <class T>
    struct value_rank { T value; int rank; };

    /// @brief Минимальное значение и номер процесса, на котором достигается
    static value_rank<double> argmin(double value) {
        return {.value=value, .rank=rank()};
    }

    /// @brief Максимальное значение и номер процесса, на котором достигается
    static value_rank<double> argmax(double value) {
        return {.value=value, .rank=rank()};
    }

    /// @}

    /// @{
    /// @name Коллективные обменные операции

    /// @brief Отправляет сообщение с "корневого" процесса всем процессам.
    /// @param root Ранг "корневого" процесса.
    /// @param value Величина шаблонного типа по ссылке, переменная также
    /// принимает пересланное значение.
    template <class T>
    static void broadcast(int root, T& value) { }

    /// @brief Коллективная операция. Собирает величину value со всех
    /// процессов, записывает в массив и раздает этот массив всем процессам.
    template <class T>
    static std::vector<T> all_gather(const T& value) { return {value}; }

    /// @brief Коллективная операция. Собирает величину value со всех процессов,
    /// записывает в массив и раздает этот массив всем процессам сети.
    template <class T>
    static void all_gather(const T& value, std::vector<T>& values) {
        values[0] = value;
    }

    /// @brief Коллективная операция. Размеры буфферов send и recv совпадают
    /// и равны числу процессов. Каждая величина из send отправляется своему
    /// процессу, каждая величина в recv получается с определенного процесса.
    template <class T>
    static void all_to_all(const std::vector<T>& send, std::vector<T>& recv) {
        recv = send;
    }

    /// @}
};

} // namespace zephyr::utils

#endif