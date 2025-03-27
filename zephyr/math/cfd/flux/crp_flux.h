#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/matter/materials.h>


namespace zephyr::math {

/// @brief Решатель композитной задачи Римана.
/// Не является наследником класса NumFlux для классических решателей,
/// поскольку функция расчета потока не является автомодельной (зависимость
/// от шага интегрирования), а также зависит от дополнительного параметра
/// delta -- отступ от грани.
class CrpFlux {
public:
    /// @brief Классическая составная задача Римана. Левая ячейка содержит
    /// вектора состояния zLA, zLB, правая ячейка только состояние zRB.
    /// Граница раздела между zLA и zLB располагается на расстоянии
    /// delta от грани между ячейками, через которую вычисляется поток.
    /// @param zLA, zLB, zRB Векторы состояния
    /// @param mixture PT-смесь, требуется для оценки скорости звука
    /// @param delta Отступ от грани ячейки влево (delta > 0)
    /// @param dt Шаг интегрирования по времени
    static mmf::Flux classic(const mmf::PState &zLA, const mmf::PState &zLB, const mmf::PState &zRB,
                             const phys::MixturePT &mixture, double delta, double dt);

    /// @brief Инвертированная составная задача Римана. Левая ячейка содержит
    /// единственный вектор состояния zLA, правая ячейка состояния zRA и zRB.
    /// Граница раздела между zRA и zRB располагается на расстоянии
    /// delta от грани между ячейками, через которую вычисляется поток.
    /// @param zLA, zLA, zRB Вектора состояния
    /// @param mixture PT-смесь, требуется для оценки скорости звука
    /// @param delta Отступ от грани ячейки вправо (delta > 0)
    /// @param dt Шаг интегрирования по времени
    static mmf::Flux inverse(const mmf::PState &zLA, const mmf::PState &zRA, const mmf::PState &zRB,
                             const phys::MixturePT &mixture, double delta, double dt);
};

} // namespace zephyr::math