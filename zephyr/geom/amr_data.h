#pragma once

namespace zephyr { namespace geom {

/// @brief Информация о ячейке адаптивной сетки
struct AmrData {
    int base_id;  ///< Индекс среди базовых ячеек
    int z;        ///< Индекс ячейки на z-кривой
    int next;     ///< Индекс новой ячейки (в алгоритмах)

    short level;  ///< Уровень адаптации (0 для базовой)
    short flag;   ///< Желаемый флаг адаптации

	/// @brief Конструктор по-умолчанию
	AmrData() : base_id(0), z(0), next(-1), level(0), flag(0) { }
};

} // geom
} // zephyr
