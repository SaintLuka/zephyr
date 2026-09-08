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

    Variables variables;  ///< Список переменных на запись
    VtuOptions options;   ///< Опции записи

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
    void open(std::string_view filename);
    void open(std::string_view filename, bool distributed);
    void open(std::string_view filename, std::string_view directory);
    void open(std::string_view filename, std::string_view directory, bool distributed);

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

    bool           open_;         ///< Открыт ли PVD файл?
    bool           distributed_;  ///< Общий PVD при использовании MPI?
    std::string    filename_;     ///< Имя файла без расширения
    std::string    fullname_;     ///< Абсолютное имя файла без расширения
    std::streamoff pos_;          ///< Указатель на позицию в файле
    std::size_t    counter_;      ///< Счетчик записанных временных шагов
};

} // namespace zephyr::io