#pragma once

#include <zephyr/io/vtu_file.h>

namespace zephyr::utils { class Json; }

namespace zephyr::io {

/// @brief Запись серии VTU-файлов и соответствующего PVD-файла.
/// @details Данный класс позволяет сохранять файлы как в однопроцессорном
/// режиме, так и при параллельном счете.
class PvdFile {
public:
    // Переменные класса имеют публичный доступ, сохранение можно выполнить
    // после настройки всех параметров

    Variables variables;   ///< Список переменных на запись

    /// @brief Для адаптивных сеток: записывать ячейки как полигоны? В обратном
    /// случае ячейки пишутся как простые четырехугольники/шестигранники,
    /// то есть с "висящими" узлами.
    ///
    /// Для трёхмерных сеток: интерпретировать ячейки как многогранники общего вида.
    /// Сохраняет отдельно грани. Работает долго, включать при необходимости.
    /// Необходимо использовать при записи сетки с многогранниками.
    bool polyhedral = false;

    /// @brief Сохранять уникальные вершины, актуально для EuMesh. Если сетка
    /// часто перестраивается, то функция выполняется достаточно долго.
    bool unique_nodes = false;

    /// @brief Пустой конструктор, не создает PVD файл, после создания
    /// экземпляра класса требуется вызов функции PvdFile::open.
    PvdFile();

    /// @brief Открывает PVD файл для записи, записывает заголовок.
    /// @details Заполняет приватные поля класса. Вызов функции требуется,
    /// если использовался тривиальный конструктор PvdFile().
    /// @param filename Короткое имя файла, без расширения.
    /// @param directory Абсолютный или относительный путь. Можно задать
    /// пустую строку, по умолчанию directory = "output".
    /// @param distributed Один PVD на несколько процессов? При использовании
    /// MPI по умолчанию ставится distributed = true.
    void open(const char* filename);
    void open(const char* filename, bool distributed);
    void open(const char* filename, const char* directory);
    void open(const char* filename, const char* directory, bool distributed);

    void open(const std::string& filename);
    void open(const std::string& filename, bool distributed);
    void open(const std::string& filename, const std::string& directory);
    void open(const std::string& filename, const std::string& directory, bool distributed);

    /// @brief Открывает PVD файл для записи, записывает заголовок.
    /// @details Заполняет приватные поля класса. Вызов функции требуется,
    /// если использовался конструктор PvdFile без указания имени файла.
    /// Вызывает тривиальный конструктор PvdFile(), а затем одну из функций
    /// PvdFile::open(...), смотри аргументы функций PvdFile::open(...).
    template <class... Args>
    explicit PvdFile(Args&&... args) : PvdFile() {
        open(std::forward<Args>(args)...);
    }

    /// @brief Создание файла из .json конфигурации
    explicit PvdFile(const utils::Json& config);

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(mesh::EuMesh& mesh, double timestep);

    /// @brief Записать хранилище (или часть, при распределенном счете) в один
    /// файл VTU (или набор VTU), затем обновить PVD файл.
    /// Используется функция VtuFile::write
    void save(mesh::AmrCells& elements, double timestep);

private:
    std::string get_filename() const;

    void update_pvd(double timestep);

    bool           m_open;         ///< Открыт ли PVD файл?
    bool           m_distributed;  ///< Общий PVD при использовании MPI?
    std::string    m_filename;     ///< Имя файла без расширения
    std::string    m_fullname;     ///< Абсолютное имя файла без расширения
    std::streamoff m_pos;          ///< Указатель на позицию в файле
    std::size_t    m_counter;      ///< Счетчик записанных временных шагов
};

} // namespace zephyr::io