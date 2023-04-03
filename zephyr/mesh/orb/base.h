#pragma once
#include <zephyr/data/storage.h>
#include <zephyr/geometry/domain.h>
#include <zephyr/network/decomposition/type/list.h>
#include <zephyr/network/mpi/network.h>
#include <zephyr/performance/load/storage-size.h>
#ifdef ZEPHYR_ENABLE_YAML
#include <yaml-cpp/yaml.h>
#endif
#ifdef ZEPHYR_ENABLE_MULTITHREADING
#include <zephyr/multithreading/thread-pool.h>
#endif

#include <zephyr/performance/timer/stopwatch.h>

namespace zephyr { namespace network { namespace decomposition {

using ::zephyr::data::Storage;
using ::zephyr::data::Vector3d;
using ::zephyr::geometry::Domain;
using ::zephyr::data::item;
using Measurer = ::zephyr::performance::load::Base;
using ::zephyr::performance::load::StorageSize;
using ::zephyr::network::mpi::Network;
#ifdef ZEPHYR_ENABLE_MULTITHREADING
using ::zephyr::multithreading::ThreadPool;
using ::zephyr::multithreading::dummy_pool;
#endif

using ::zephyr::data::coords;
using ::zephyr::data::neibsSearchRadius;

using ::zephyr::performance::timer::Stopwatch;

using load_t = double;

static std::vector<long long int> dummy_index;

/** 
	\brief 
		\~russian Базовый класс декомпозиции данным по MPI-процессам.
		\~english Base data decomposition by MPI processes.
		\~
*/
class Base {
public:
	/** 
		\brief 
			\~russian Конструктор класса.
			\~english Class constructor.
			\~
		\param[in] net
			\~russian сеть MPI-процессов.
			\~english MPI network.
			\~
		\param[in] elements
			\~russian хранилище с элементами, для которых
				осуществляется декомпозиция.
			\~english Storage with elements for which 
				decomposition is performed.
			\~
		\param[in] domain
			\~russian геометрия области, на которой осуществляется 
				декомпозиция.
			\~english geometry domain for decomposition 
				to be performed on.
			\~
		\param[in] measurer
			\~russian измеритель нагрузки.
			\~english load measurer.
			\~
	*/
	Base(
		Network&  net,
		Storage&  elements,
		Domain&   domain,
		Measurer& measurer
        if_multithreading(, ThreadPool& pool)
	);

#ifdef ZEPHYR_ENABLE_YAML
    /**
        \brief
            \~russian Создает экзепляр класса-потомка.
            \~
        \param[in] config
            \~russian Конфигурация декомпозиции.
            \~english Config of decomposition.
            \~
    */
	static std::unique_ptr<Base> create(
            Network&  net,
            Storage&  elements,
            Domain&   domain,
            Measurer& measurer,
            if_multithreading(ThreadPool& pool,)
            const YAML::Node& config);
#endif

	/** \brief
		\~russian Виртуальный деструктор.
		\~english Virtual destructor.
		\~
	*/
	virtual ~Base() = default;

	/** \brief
		\~russian Обновить декомпозицию.
		\~english Update decomposition according to 
			input Storage changes.
		\~
	*/
	void update();

	/** \brief
		\~russian Обменяться элементами у границ подобластей декомпозиции.
			Элементы хранилища не меняются.
		\~english Exchange elements near edges of decomposition 
			subdomains, input Storage is not changed
		\~
		\see elements.
	*/
	void exchange();
	/** \brief
		\~russian Начать асинхронный обмен элементами у границ 
			подобластей декомпозиции.
		\~english Start asynchronic exchange of elements near edges.
		\~
		\see exchange.
	*/
	void exchange_start();
	/** \brief
		\~russian Закончить асинхронный обмен элементами у границ 
			подобластей декомпозиции.
		\~english Finish asynchronic exchange of elements near edges.
		\~
		\see exchange.
	*/
	void exchange_end();
	/** \brief
		\~russian Перевести элементы в глобальную систему координат.
		\~english Put elements into global coordinate system.
		\~
	*/
	void to_global_CS();
	/** \brief
		\~russian Перевести элементы в локальную систему координат.
		\~english Put elements into local coordinate system.
		\~
	*/
	void to_local_CS();

	/** \brief
		\~russian Ссылка на элементы, хранящиеся в данном процессе.
		\~english Reference to elements that are stored by the process.
		\~
	*/
	Storage& inner_elements();
	/** \brief
		\~russian Ссылка на элементы, которые не хранятся 
			в данном процессе, но используются им.
		\~english Reference to elements that are not stored 
			by the process, but are used by it.
		\~
	*/
	Storage& outer_elements();
	
	/** \brief
		\~russian Функция запуска перестановки элементов.
		\~english Function to call elements permutation function.
		\~
	*/
	void signal_to_permute(std::vector<std::size_t>& perm);

	/** \brief
		\~russian Собрать все элементы хранилища в одном MPI-процессе.
		\~english Gather all Storage elements in a single MPI process.
		\~
		\param[in] rank
		\~russian ранг MPI-процесса.
		\~english MPI-process rank.
		\~
	*/
	void collect_at(std::size_t rank);
	/** \brief
		\~russian Назначить измеритель нагрузки.
		\~english Set load measurer.
		\~
		\param[in] lm
		\~russian измеритель нагрузки.
		\~english load measurer.
		\~
	*/
	void set_load_measurer(Measurer& lm);

	/** \brief
		\~russian Получить измеритель нагрузки.
		\~english Get load measurer.
		\~
	*/
	Measurer& load_measurer() const;
	/** \brief
		\~russian Получить сеть MPI-процессов.
		\~english Get MPI network.
		\~
	*/
	Network& network() const;
	/** \brief
		\~russian Получить средние значения по списку полей величин.
		\~english Evaulate average fields.
		\~
		\param[in] fields
		\~russian список полей Storage в виде массива строк.
		\~english Storage fields as strings.
		\~
		\param[in] index
		\~russian индексы элементов для усреднения.
		\~english element indices to account.
		\~
	*/
	Storage average(
		std::vector<std::string> fields,
		std::vector<long long int>& indices = dummy_index
	);
	/** \brief
		\~russian Получить средние значения по списку полей величин.
		\~english Evaulate average fields.
		\~
		\param[in] fields
		\~russian список полей Storage в виде массива строк.
		\~english Storage fields as strings.
		\~
		\param[in] index
		\~russian указатель на массив индексов элементов для усреднения.
		\~english pointer to array of element indices to account.
		\~
		\param[in] size
		\~russian число элементов массива index.
		\~english number of elements in index array.
		\~
	*/
	Storage average(
		std::vector<std::string> fields,
		std::pair<long long int*, std::size_t> indices = {nullptr, 0}
	);

    /** \~russian Получить "центр" части декомпозиции.
        \~
    */
    virtual Vector3d center() = 0;

    /// Создать списки aliens для отправки
    /// Формирует список tourists_labels_out
    void prepare_aliens();

    /// Всякие таймеры
    Stopwatch& prep_timer()  { return m_prep_timer; }
    Stopwatch& send_timer()  { return m_send_timer; }
    Stopwatch& send1_timer() { return m_send1_timer; }
    Stopwatch& send2_timer() { return m_send2_timer; }
    Stopwatch& recv_timer()  { return m_recv_timer; }
    Stopwatch& wait_timer()  { return m_wait_timer; }

protected:

    /// Улучшенный интерфейс (в разработке)

    virtual void collect_local_info() = 0;

    virtual void collect_aliens_info() = 0;

    /// @brief Важная функция. Определить ранг процесса, которому принадлежит
    /// произвольная точка v. На практике точка v обычно является положением
    /// частицы или центром расчетной ячейки.
    /// Функция используется для перераспределения элементов между процессами.
    /// @param v Произвольная точка всей области
    virtual int rank(const Vector3d& v) const = 0;

    /// @brief Важная функция. Определить, принадлежит ли некоторая точка v
    /// данного (!) процесса в некоторой окрестности процесса с рангом
    /// neib_rank. Функция используется для формирования списков aliens.
    /// В качестве радиуса поиска обычно выбирается максимальный радиус
    /// поиска для элементов с данного процесса и с соседнего (neib_rank)
    /// @param v Точка из подобласти данного процесса
    /// @param neib_rank Ранг процесса, в окрестности которого ищется точка
    /// @return true если точка принадлежит окрестности
    virtual bool is_near(const Vector3d& v, int neib_rank) const = 0;


    /** \~russian Пометить элемент для обмена с соседними процессами.
        \~english Mark element for exchange with neighboring processes.
        \~
        \param[in] it
            \~russian элемент.
            \~english element.
            \~
        \param[in] begin
            \~russian начало хранилища элементов.
            \~english elements Storage begin.
            \~
        \see tourists_labels_out.
    */
    void mark_tourists_out(Storage::Item it, Storage::Item begin);

    /// Переслать aliens
    //void exchange_aliens();

    /// Собрать списки migrants
    /// Сформировать migrants_out и to_remove
    void prepare_migrants();

    virtual void balancing(const std::vector<double>& w) = 0;

    /** \~russian Отдать элементы, которые больше не принадлежат
            этому процессу, другим процессам.
        \~english Give away elements which no longer belong to the diagram
            subdomain to other processes.
        \~
    */
    void exchange_migrants();

    /// Обменяться migrants
    // void exchange_migrants();

    /// Перестроить диаграммы или что там есть у конкретного метода,
    /// без изменения каких-либо списков мигрантов и туристов.
    //virtual void balancing() = 0;

    /// Начальная штука
    //void redistribute();

    // Отправляет данные без изменений списков
    // exchange_start();

    // Получает данные без каких-либо изменений списков
    // exchange_end();

    // Сейчас делает лишний prepare, неочевидное поведение,
    // должно быть просто равноценно последовательным вызовам
    // exchange_start() и exchange_end()
    // exchange();

    void redistribute();

    // В некотором смысле функция перегружена, очень сложное
    // поведение: изменение декомпозиции, поиск и маркировка
    // приграничных элементов, затем пересылка элементов.
    // update() {
    //    balancing() // virtual
    //
    //    redistribute() //
    //}

    /** \~russian Инициализировать список соседей процесса.
        \~english Initialize neighbors array.
        \~
    */
    void init_neibs();

	/** \brief
		\~russian Инициализировать переменные, связанные с
			обменным слоем.
		\~english Initialize variables related to an exchange layer.
		\~
	*/
	void init_tourists();
	/** \brief
		\~russian Инициализировать переменные, связанные с
			передачей элементов другим процессам.
		\~english Initialize variables related to transferring
			elements.
		\~
	*/
	void init_migrants();
	/** \brief
		\~russian Осуществить перестановку элементов в массиве
			tourists_labels_out на основании входного вектора.
		\~english Perform permutation of the elements in 
			tourists_labels_out according to the input Vector3d.
		\~
		\see tourists_labels_out.
	*/
	void permute_tourists_labels_out(std::vector<std::size_t>& perm);

    /** \~russian Поместить в период элементы из хранилища.
        \~english Fit elements in period of domain.
        \~
        \param[in] elements
        \~russian хранилище элементов.
        \~english Storage with elements.
        \~
    */
    void fit_in_period_for(Storage& elements);

    /** \~russian Разместить элементы ближе вокруг центра диаграммы,
            если позволяет периодичность области.
        \~english Place elements nearer to the center of the diagram
             if periodicity of the domain allows.
        \~
        \param[in] elements
        \~russian хранилище элементов.
        \~english Storage with elements.
        \~
    */
    void remove_period_for(Storage& elements);

    /** \~russian Поместить элемент в период.
        \~english Fit an element in period.
        \~
        \param[in] i
            \~russian итератор элемента.
            \~english an iterator to an element.
            \~
    */
    void fit_in_period(Storage::Item i);

    /** \~russian Разместить элемент ближе к центру диаграммы,
            если позволяет периодичность области.
        \~english Place an element nearer to the center of the diagram
             if periodicity of the domain allows.
        \~
        \param[in] element
            \~russian итератор элемента.
            \~english an iterator to an element.
            \~
    */
    void remove_period(Storage::Item element);

	Network& net;             ///< Сеть MPI-процессов
	Storage& m_locals;       ///< Локальное хранилище данных
    Domain m_domain;          ///< Геометрия расчетной области
    Measurer* m_measurer;     ///< Измеритель нагрузки
#ifdef ZEPHYR_ENABLE_MULTITHREADING
    ThreadPool& m_pool;       ///< Ссылка на пул тредов для многопоточности
#endif

    int m_rank;               ///< Ранг текущего процесса

	bool has_exchange_layer;  ///< Есть ли обменный слой

    /// @brief Список соседей текущего процесса, то есть ранги процессов,
    /// с которыми происходит обмен данными приграничными элементами
    /// @details Список соседей так же требуется в алгоритмах балансировки
    /// с использованием диаграммы Вороного.
    std::vector<int> m_neibs;

    /// @brief Массив индексов элементов из обменного слоя, которые необходимо
    /// передавать другим процессам.
    std::vector<std::vector<size_t>> m_tourists_labels;

    /// @brief Хранище для элементов из обменного слоя данного процесса.
	Storage m_tourists;

    /// @brief Хранилище элементов из обменного слоя, которые переданы
    /// от других процессов.
	Storage m_aliens;

	std::vector<size_t>      m_send_count;    ///< Число ячеек для отправки
    std::vector<MPI_Request> m_send_request;  ///< MPI-запрос на отправку
    std::vector<MPI_Status>  m_send_status;   ///< MPI-статус отправки

    std::vector<size_t>      m_recv_count;    ///< Число ячеек для отправки
    std::vector<MPI_Request> m_recv_request;  ///< MPI-запрос на получение
    std::vector<MPI_Status>  m_recv_status;   ///< MPI-статус получения




    /// @brief Хранище для элементов, переданных от других процессов
    std::vector<Storage> migrants_in;

    /// @brief Хранище для элементов, которые должны быть переданы другим
    /// процессам
    std::vector<Storage> migrants_out;

    /// @brief Список элементов, которые будут удалены.
    std::vector<std::size_t> to_remove;


    Stopwatch m_prep_timer;
	Stopwatch m_send_timer;
    Stopwatch m_send1_timer;
    Stopwatch m_send2_timer;
	Stopwatch m_recv_timer;
	Stopwatch m_wait_timer;
};


} /*decomposition*/ } /*network*/ } /*zephyr*/