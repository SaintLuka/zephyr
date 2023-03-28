#pragma once

#if IM_ENABLE_MPI

#include <mpi.h>
#include <control/mpi_wrapper.h>

#include <core/layer/layer_control.h>

/**
 * Класс для обмена массивами произвольных структур между случайными процессами.
 */
class MpiCommutator {
public:

    /**
     * @brief Конструктор класса.
     *
     * @param need_to_send Вектор содержит информацию о том,
     * какому процессу необходимо передавать данные.
     * need_to_send[r] - необходимо передать данные процессу r?
     */
    explicit MpiCommutator(const vector<bool>& need_to_send) {
        vector<int> send = vector<int>(mpi::usize(), 0);
        for (size_t r = 0; r < need_to_send.size(); ++r) {
            send[r] = need_to_send[r] ? 1 : 0;
        }

        vector<int> recv = vector<int>(mpi::usize(), 0);
        for (int r = 0; r < mpi::size(); ++r) {
            MPI_Gather(&send[r], 1, MPI_INT, recv.data(), 1, MPI_INT, r, MPI_COMM_WORLD);
        }

        m_need_to_send = vector<bool>(mpi::usize(), false);
        m_need_to_recv = vector<bool>(mpi::usize(), false);
        for (int r = 0; r < mpi::size(); ++r) {
            m_need_to_send[r] = send[r] > 0;
            m_need_to_recv[r] = recv[r] > 0;
        }
        m_need_to_send[mpi::rank()] = false;
        m_need_to_recv[mpi::rank()] = false;

        // dbg_info();
    }


    /**
     * @brief Обменивается данными между процессами.
     * Все полученные данные собираются в единый буффер.
     *
     * @tparam T Тип передаваемых данных. Структура/класс с простой сериализацией
     * (любой тип приводится к MPI_BYTE)
     *
     * @param s_buffers Набор массивов для отпавки каждому из процессов.
     * Массив s_buffers[r] будет отправлен r-ому процессу.
     *
     * @return Единый массив, в котором собраны все полученные данные
     */
    template<class T>
    vector<T> send_to_common_buffer(const vector<vector<T>> &s_buffers) const {
        vector<int> s_sizes;
        vector<int> r_sizes;
        calc_sizes(s_buffers, s_sizes, r_sizes);

        vector<int> offsets(mpi::usize() + 1, 0);
        for(int r = 0; r < mpi::size(); ++r) {
            offsets[r + 1] = offsets[r] + r_sizes[r];
        }

        vector<T> r_buffer(static_cast<unsigned int>(offsets.back()));

        vector<MPI_Request> recv_request(mpi::usize());
        vector<MPI_Request> send_request(mpi::usize());

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Irecv(&r_buffer[offsets[r]], r_sizes[r] * sizeof(T), MPI_BYTE, r,
                          TAGS::COUNT_TAG, MPI_COMM_WORLD, &recv_request[r]);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Isend(s_buffers[r].data(), s_sizes[r] * sizeof(T), MPI_BYTE, r,
                          TAGS::COUNT_TAG, MPI_COMM_WORLD, &send_request[r]);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Wait(&recv_request[r], MPI_STATUS_IGNORE);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Wait(&send_request[r], MPI_STATUS_IGNORE);
            }
        }

        return r_buffer;
    }


    /**
     * @brief Обменивается данными между процессами.
     * Данные, полученные от каждого процесса, хранятся в различных массивах.
     *
     * @tparam T Тип передаваемых данных. Структура/класс с простой сериализацией
     * (любой тип приводится к MPI_BYTE)
     *
     * @param s_buffers Набор массивов для отпавки каждому из процессов.
     * Массив s_buffers[r] будет отправлен r-ому процессу.
     *
     * @return Набор векторов, которые содержат полученные данные.
     * out[r] - массив, полученный от r-ого процесса.
     */
    template<class T>
    vector<vector<T>> send_to_separated_buffers(const vector<vector<T>> &s_buffers) const {

        vector<int> s_sizes;
        vector<int> r_sizes;
        calc_sizes(s_buffers, s_sizes, r_sizes);

        vector<vector<T>> r_buffers(mpi::usize());
        for(int r = 0; r < mpi::size(); ++r) {
            r_buffers[r] = vector<T>(static_cast<unsigned int>(r_sizes[r]));
        }

        vector<MPI_Request> recv_request(mpi::usize());
        vector<MPI_Request> send_request(mpi::usize());

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Irecv(r_buffers[r].data(), r_sizes[r] * sizeof(T), MPI_BYTE, r,
                          TAGS::COUNT_TAG, MPI_COMM_WORLD, &recv_request[r]);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Isend(s_buffers[r].data(), s_sizes[r] * sizeof(T), MPI_BYTE, r,
                          TAGS::COUNT_TAG, MPI_COMM_WORLD, &send_request[r]);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Wait(&recv_request[r], MPI_STATUS_IGNORE);
            }
        }

        for(int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Wait(&send_request[r], MPI_STATUS_IGNORE);
            }
        }

        return r_buffers;
    }


private:

    /**
     * @brief Определяет размеры отправляемых/получаемых данных
     *
     * @tparam T Тип передаваемых данных
     *
     * @param s_buffers Набор массивов для отпавки каждому из процессов.
     * Массив s_buffers[r] будет отправлен r-ому процессу.
     *
     * @param s_sizes Выходной массив.
     * После выполнения содержит размеры массивов,
     * отправляемых каждому из процессов.
     *
     * @param r_sizes Выходной массив.
     * После выполнения содержит размеры массивов,
     * получаемых от каждого из процессов.
     */
    template<class T>
    void calc_sizes(const vector<vector<T>>& s_buffers, vector<int>& s_sizes, vector<int>& r_sizes) const {

        s_sizes.resize(mpi::usize());
        for (int r = 0; r < mpi::size(); ++r) {
            s_sizes[r] = static_cast<int>(s_buffers[r].size());
        }

        r_sizes.resize(mpi::usize());

        vector<MPI_Request> count_recv_request(mpi::usize());
        vector<MPI_Request> count_send_request(mpi::usize());

        for (int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Irecv(&r_sizes[r], 1, MPI_INT, r, TAGS::COUNT_TAG, MPI_COMM_WORLD, &count_recv_request[r]);
            }
        }

        for (int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Isend(&s_sizes[r], 1, MPI_INT, r, TAGS::COUNT_TAG, MPI_COMM_WORLD, &count_send_request[r]);
            }
        }

        for (int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_recv[r]) {
                MPI_Wait(&count_recv_request[r], MPI_STATUS_IGNORE);
            }
        }

        for (int r = 0; r < mpi::size(); ++r) {
            if (m_need_to_send[r]) {
                MPI_Wait(&count_send_request[r], MPI_STATUS_IGNORE);
            }
        }
    }


    /**
     * Выводит информацию о передачах
     */
    void dbg_info() {
        for (int r = 0; r < mpi::size(); ++r) {
            mpi::barrier();
            if (mpi::rank() == r) {
                std::cout << "  hi from rank " << mpi::rank() << "\n";
                std::cout << "    send to: ";
                for (int r2 = 0; r2 < mpi::size(); ++r2) {
                    if (m_need_to_send[r2])
                        std::cout << r2 << ", ";
                }
                std::cout << "\n";

                std::cout << "    wait from: ";
                for (int r2 = 0; r2 < mpi::size(); ++r2) {
                    if (m_need_to_recv[r2] > 0)
                        std::cout << r2 << ", ";
                }
                std::cout << "\n";
            }
        }
    }


    /** Каким процессам отправляются данные */
    vector<bool> m_need_to_send;


    /** От каких процессов получаем данные */
    vector<bool> m_need_to_recv;
};

#endif
