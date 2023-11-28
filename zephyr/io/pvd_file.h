#pragma once

#include <zephyr/configuration.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr::io {

/// @class Класс для записи серии VTU файлов и соответствующего PVD файла.
/// @details Данный класс позволяет сохранять файлы как в однопоточном режиме,
/// так и при параллельном счете.
class PvdFile {
public:
    /// @brief Набор public переменных: фильр, список переменных и тип
    /// сохранения сетки могут настраиваться пользователем класса.

    Variables     variables; ///< Список переменных на запись
    ComplexFilter filter;    ///< Функтор, возвращает true для нужных ячеек
    bool          hex_only;  ///< Записывать как четырехугольники/шестигранники

    /// @brief Пустой конструктор, не создает PVD файл, после создания
    /// экземпляра класса требуется вызов функции PvdFile::open.
    PvdFile();

    /// @brief Конструктор класса. При создании экземпляра класса сразу
    /// открывается PVD файл и записываететя заголовок.
    /// @param filename Короткое имя файла, без расширения
    /// @param directory Абсолютный или относительный путь
    explicit PvdFile(const std::string& filename, const std::string& directory = "");

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    /// @brief Пустой конструктор, не создает PVD файл, после создания
    /// экземпляра класса требуется вызов функции PvdFile::open.
    explicit PvdFile(Network& network);

    /// @brief Конструктор класса
    /// @param network Необходимо передать сеть, если требуется сохранять набор файлов.
    /// @param filename Короткое имя файла, без расширения
    /// @param directory Абсолютный или относительный путь
    PvdFile(Network& network, const std::string& filename, const std::string& directory = "");
#endif

    /// @brief Открывает PVD файл для записи, записывает заголовок.
    /// @details Заполняет приватные поля класса. Вызов функции требуется, если
    /// использовался конструктор PvdFile без указания имени файла.
    /// @param filename Короткое имя файла, без расширения
    /// @param directory Абсолютный или относительный путь
    void open(const std::string& filename, const std::string& directory = "");

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(mesh::EuMesh& mesh, double timestep);

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(AmrStorage& elements, double timestep);

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(mesh::LaMesh& mesh, double timestep);

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(CellStorage& cells, NodeStorage& nodes, double timestep);

private:

    std::string get_filename() const;

    void update_pvd(double timestep);

    bool           m_open;     ///< Открыт ли pvd файл
    std::string    m_filename; ///< Имя файла без расширения
    std::string    m_fullname; ///< Абсолютное имя файла без расширения
    std::streamoff m_pos;      ///< Указатель на позицию в файле
    std::size_t    m_counter;  ///< Счетчик записанных временных шагов

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    Network& net;
#endif
};

} // namespace zephyr::io