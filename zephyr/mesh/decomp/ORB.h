#pragma once

#include <zephyr/mesh/decomp/orb/blocks.h>
#include <zephyr/mesh/decomp/decomposition.h>

namespace zephyr::mesh::decomp {

class ORB : public Decomposition {
public:
    /// @brief Умный указатель на экземпляр класса
    using Ptr = std::shared_ptr<ORB>;
    using Ref = const std::shared_ptr<ORB>&;

    /// @brief Параметры декомпозиции
    struct params {
        bool   newton = false;  ///< Использовать метод Ньютона
        double mobility = 0.1;  ///< Скрость смещения генераторов (0, 0.5)
    };

    /// @brief Декмопозиция для прямоугольной области
    /// @param domain Прямоугольная/кубическая область
    /// @param type Тип декмопозиции: X, XY, YX, XYZ, ...
    /// @param size Количество блоков
    /// @param p Опции балансировки
    ORB(Box domain, const std::string& type, int size,
        const params& p = {.newton = false, .mobility = 0.1});

    /// @brief Конструктор с автоматической декомпозицией по размерам области
    /// @param box Домен
    /// @param type Тип декомпозиции
    /// @param size Число блоков
    /// @param nx Число блоков по первой координате
    ORB(Box domain, const std::string& type, int size,
        int nx, const params& p = {.newton = false, .mobility = 0.1});

    /// @brief Конструктор с автоматической декомпозицией по размерам области
    /// @param box Домен
    /// @param type Тип декомпозиции
    /// @param size Число блоков
    /// @param ny Число блоков по второй координате, для двумерной декомпозиции
    /// это полное описание декомпозиции
    ORB(Box domain, const std::string& type, int size, const std::vector<int>& ny,
        const params& p = {.newton = false, .mobility = 0.1});

    /// @brief Создание указателя
    static ORB::Ptr create(const Box& domain, const std::string& type, int size,
                           const params& p = {.newton = false, .mobility = 0.1}){
        return std::make_shared<ORB>(domain, type, size, p);
    }

    /// @brief Создание указателя
    static ORB::Ptr create(const Box& domain, const std::string& type, int size, int nx,
                           const params& p = {.newton = false, .mobility = 0.1}) {
        return std::make_shared<ORB>(domain, type, size, nx, p);
    }

    /// @brief Создание указателя
    static ORB::Ptr create(const Box& domain, const std::string& type, int size,
                           const std::vector<int>& ny,
                           const params& p = {.newton = false, .mobility = 0.1}){
        return std::make_shared<ORB>(domain, type, size, ny, p);
    }

    /// @brief Определить ранг процесса, которому принадлежит точка v.
    /// На практике точка v обычно является положением центра ячейки
    int rank(const Vector3d& v) const;

    int rank(AmrStorage::Item& elem) const override;

    // set функции

    void use_newton(bool val = true);

    void set_mobility(double val);

    // get функции

    bool newton() const;

    double mobility() const;

    void balancing(const std::vector<double>& w) final;

private:

    /// @brief Использовать метод Ньютона
    bool m_newton;

    /// @brief Скорость сдивга (0, 0.5)
    double m_mobility;

    /// @brief Блочная структура
    Blocks m_blocks;
};

} // namespace zephyr::mesh::decomp
