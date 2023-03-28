#pragma once
#include <vector>
#include <zephyr/network/mpi/node.h>

namespace zephyr { namespace network { namespace mpi {
/** 
	\brief 
		\~russian Класс абстрактной сети вычислительных узлов, 
		          внутри которой можно пересылать данные.
		\~english An abstract class for computational nodes network,
		          among which data can be exchanged.
		\~
*/
class Network {
public:
	std::VECTOR<Node> node; /**<
		\~russian Список вычислительных узлов для использования.
		\~english A list of computational nodes for usage.
		\~
	*/
private:
	MPI_Comm _comm; /**<
		\~russian Коммуникатор MPI-процессов данной сети.
		\~english An MPI communicator of this Network.
		\~
	*/
	MPI_Group _group; /**<
		\~russian Группа MPI-процессов данной сети.
		\~english An MPI group of this Network.
		\~
	*/
    int _rank; /**<
		\~russian Номер MPI-процесса, которому соответствует
		          вычислительный узел.
		\~english An MPI-process rank attached to Node.
		\~
	*/
	int _size; /**<
		\~russian Число MPI-процессов в данной сети.
		\~english A number of MPI processes in this Network.
		\~
	*/
	std::VECTOR<Network*> _subnet; /**<
		\~russian Список подсетей данной сети.
		\~english A list of subnetworks of this Network.
		\~
	*/
	Node dummy; /**<
		\~russian Узел-заглушка для MPI-процесса, который 
		          находится вне данной сети.
		\~english A dummy node for an MPI process outside
		          this Network.
		\~
	*/
public:
	/** 
		\brief 
			\~russian Запрет на использование конструктора по-умолчанию.
			\~english Prohibit using default constructor.
			\~
	*/
	Network() = delete;
	/** 
		\brief 
			\~russian Конструктор сети с помощью MPI-коммуникатора.
			\~english Network constructor using MPI-communicator.
			\~
	*/
	Network(MPI_Comm c);
	/** 
		\brief 
			\~russian Запрет на использование оператора присваивания сети.
			\~english Prohibit using assignment operator.
			\~
	*/
	Network& operator=(const Network& N) = delete;
	/** 
		\brief 
			\~russian Создать подсеть из списка порядковых номеров узлов.
			\~english Create subnetwork from a list of node numbers.
			\~
	*/
	Network& subnet(std::VECTOR<int>& ranks);
	/** 
		\brief 
			\~russian Получить ссылку на узел, где исполняется
			          данный MPI-процесс.
			\~english Get a node reference to a node where this 
			          MPI-process is executed.
			\~
	*/
	Node& self();
	/** 
		\brief 
			\~russian Получить число узлов в сети.
			\~english Get a number of nodes in network.
			\~
	*/
	std::size_t size() const;
    /**
        \brief
            \~russian Получить номер MPI-процесса, соответствующий
                      вычислительному узлу.
            \~english Get MPI-process rank related to Node.
            \~
    */
    int rank() const;
	/** 
		\brief 
			\~russian Получить MPI-коммуникатор сети.
			\~english Get an MPI communicator of network.
			\~
	*/
	MPI_Comm comm() const;
	/** 
		\brief 
			\~russian Получить MPI-группу данной сети.
			\~english Get an MPI group of this network.
			\~
	*/
	MPI_Group group() const;
	/** 
		\brief 
			\~russian Деструктор класса.
			\~english Class destructor.
			\~
	*/
	~Network();
    /**
        \brief
            \~russian Возвращает true для процесса нулевого ранга.
            \~english Return true for zero rank process.
            \~
    */
    bool is_master() const;
    /**
        \brief
            \~russian Блокирует все процессы данной сети, пока они не
                      достигнут данного вызова.
            \~english Blocks until all processes in the commutator have
                      reached this routine.
            \~
    */
	void barrier() const;
    /**
        \brief
            \~russian Коллективная операция. Минимум величины среди всех
                      процессов сети.
            \~english Collective function. Minimum of value among all processes
                      in the communicator.
            \~
    */
    template <class T>
	T min(const T& value) const {
        T glob_value;
        MPI_Allreduce(&value, &glob_value, 1, type::of<T>(), MPI_MIN, _comm);
        return glob_value;
	}
    /**
        \brief
            \~russian Коллективная операция. Покомпонентный минимум для каждого
                      элемента вектора серди всех процессов сети.
            \~english Collective function. Minimum of values among all processes
                      in the communicator.
            \~
        \param values
            \~russian Вектор double, размеры векторов должны совпадать на каждом
                      процессе сети.
        \return
            \~russian Выходной вектор имеет размер вектора в аргументе.
    */
    template <class T>
    std::VECTOR<T> min(const std::VECTOR<T>& values) const {
        std::VECTOR<T> glob_values(values.size());
        MPI_Allreduce(values.data(), glob_values.data(), int(values.size()),
                      type::of<T>(), MPI_MIN, _comm);
        return glob_values;
    }
    /**
        \brief
            \~russian Коллективная операция. Максимум величины среди всех
                      процессов сети.
            \~english Collective function. Maximum of values among all processes
                      in the communicator.
            \~
    */
    template <class T>
    T max(const T& value) const {
        T glob_value;
        MPI_Allreduce(&value, &glob_value, 1, type::of<T>(), MPI_MAX, _comm);
        return glob_value;
    }
    /**
        \brief
            \~russian Коллективная операция. Максимум величины среди всех
                      процессов сети.
            \~english Collective function. Minimum of values among all processes
                      in the communicator.
            \~
    */
    template <class T>
    T max(const std::VECTOR<T>& values) const {
        std::VECTOR<T> glob_values(values.size());
        MPI_Allreduce(values.data(), glob_values.data(), int(values.size()),
                      type::of<T>(), MPI_MAX, _comm);
        return glob_values;
    }
    /**
        \brief
            \~russian Коллективная операция. Сумма величин со всех процессов сети.
            \~english Collective function. Sum of values from all processes
                      in the communicator.
            \~
    */
    template <class T>
    T sum(const T& value) const {
        T glob_value;
        MPI_Allreduce(&value, &glob_value, 1, type::of<T>(), MPI_SUM, _comm);
        return glob_value;
    }
    /**
        \brief
            \~russian Коллективная операция. Покомпонентная сумма элементов
                      векторов со всех процессов сети.
            \~
        \param values
            \~russian Вектор double, размеры векторов должны совпадать на каждом
                      процессе сети.
            \~
        \return
            \~russian Выходной вектор имеет размер вектора в аргументе.
            \~
    */
    template <class T>
    std::VECTOR<T> sum(const std::VECTOR<T>& values) {
        std::VECTOR<T> glob_values(values.size());
        MPI_Allreduce(values.data(), glob_values.data(), int(values.size()),
                      type::of<T>(), MPI_SUM, _comm);
        return glob_values;
    }
    /**
        \brief
            \~russian Отправляет сообщение с "корневого" процесса всем
                      процессам данной сети.
            \~english Broadcasts a message from the process with rank "root"
                      to all other processes of the communicator.
            \~
        \param root
            \~russian Ранг "корневого" процесса.
            \~english Rank of the "root" process.
            \~
        \param value
            \~russian Величина шаблонного типа по ссылке, переменная также
                      принимает пересланное значение.
            \~english Value of template type by reference, it also assigns
                      the received value.
            \~
        \return
            \~russian Выходной вектор имеет размер вектора в аргументе.
            \~
    */
    template <class T>
    void broadcast(int root, T& value) const {
        MPI_Bcast((void*)&value, sizeof(value), MPI_CHAR, root, _comm);
    }
    /**
        \brief
            \~russian Коллективная операция. Собирает величину value со всех
                      процессов сети, записывает в массив и раздает этот массив
                      всем процессам сети.
            \~
    */
    template <class T>
    std::VECTOR<T> all_gather(const T& value) {
        std::VECTOR<T> values((size_t)_size);
        MPI_Allgather(&value, 1, type::of<T>(), values.data(), 1, type::of<T>(), _comm);
        return values;
    }
    /**
        \brief
            \~russian Коллективная операция. Собирает величину value со всех
                      процессов сети, записывает в массив и раздает этот массив
                      всем процессам сети.
            \~
    */
    template <class T>
    void all_gather(const T& value, std::VECTOR<T>& values) {
        values.resize((size_t)_size);
        MPI_Allgather(&value, 1, type::of<T>(), values.data(), 1, type::of<T>(), _comm);
    }
    /**
        \brief
            \~russian Коллективная операция. Размеры буфферов send и recv
                      совпадают и равны числу процессов. Каждая величина из
                      send отправляется своему процессу, каждая величина в
                      recv получается с определенного процесса.
            \~
    */
    template <class T>
    void all_to_all(const std::VECTOR<T>& send, std::VECTOR<T>& recv) {
        recv.resize((size_t) _size);
        MPI_Alltoall(
                send.data(), 1, type::of<T>(),
                recv.data(), 1, type::of<T>(),
                _comm);
    }

};
/** 
	\brief 
		\~russian Функция доступа к корневой сети MPI_COMM_WORLD.
		\~english An access function to root network MPI_COMM_WORLD.
		\~
*/
Network& main();

extern Network* _main; /**<
	\~russian Указатель на корневую сеть MPI_COMM_WORLD.
	\~english A pointer to root network MPI_COMM_WORLD.
	\~
*/
extern Network dummy; /**<
	\~russian Сеть-заглушка для использования зависимых модулей без MPI.
	\~english A dummy network to use dependent modules without MPI.
	\~
*/

}/*network*/}/*network*/}/*zephyr*/
