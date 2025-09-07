# Общие вопросы

1. Восстановить и причесать все тесты в новой SoA-парадигме.
2. Причесать код (::mesh, ::amr), обновить inline-документацию.
3. Восстановить все функции в новом коде.

# Асинхронные обмены в MPI

Сделать и протестировать скрытые пересылки. Проверить массив unique_border_indices,
а также различные варианты организации пересылок.
Попробовать впихнуть Tourism::prepare внутрь решателя?

# Гибридная параллельность

Использование тредов в MPI-алгоритмах миграции и пересылок. Параллельное
выполнение вспомогательных алгоритмов.

# AMR + MPI

Оптимизация MPI-версии адаптации (без redistribute).
Хороший тест: симметричный sm_fluid_3D на 8 MPI-процессов. Позволяет сравнить
производительность AMR с использованием чистого MPI алгоритма и чистых тредов.

## Оценка производительности

MPI-версия (8 процессов, статичная декомпозиция 2 x 2 x 2):

Refine elapsed:      138922 ms
  Balance flags:      74520 ms
    Base restrictions: 1357 ms
    Setup vicinity:    5906 ms
    Flag balancing:   67255 ms
  Apply flags:        64401 ms
    Statistics:          42 ms
    Positions:          888 ms
    Geometry:          6095 ms
    Send/Recv:         1451 ms
    Connections:       4847 ms
    Build aliens 1:   13724 ms
    Build aliens 2:   27297 ms
    Remove undefined:     10052 ms
      Create SwapList:   72 ms
      Set mapping:       34 ms
      Send/Recv next:   518 ms
      Change adjacent: 8814 ms
      Swap elements:    607 ms

Многопоточная версия (8 тредов, TBB):

Refine elapsed:       27132 ms
  Balance flags:       3654 ms
    Base restrictions: 1666 ms
    Sort by level:      941 ms
    Round 1 elapsed:    416 ms (avg retain: 40643)
    Round 2 elapsed:    273 ms (avg coarse: 36818)
    Round 3 elapsed:    173 ms
    Round 4 elapsed:    180 ms
  Apply flags:        23478 ms
    Statistics:         156 ms
    Positions:         2834 ms
    Geometry:          6540 ms
    Connections:       4827 ms
    Remove undef:      9087 ms
      Create SwapList:  223 ms
      Set mapping:       81 ms
      Change adjacent: 7898 ms
      Swap elements:    881 ms

Скорость однопроцессорной fast_balancing поражает.

То есть производительность MPI рубится преимущественно на двух пунктах:
  1. Балансировка флагов (алгоритм slow_balancing).
  2. Двойное построение обменных слоёв build_aliens.
Непосредственно MPI-пересылки влияют не так сильно.

Slow balancing заменить на fast_balancing. В первом приближении fast_balancing
не отличается от однопроцессорной версии, просто требует большого количества
пересылок. Затем попробовать убрать часть пересылок, которые не влияют на
выполнение алгоритма. Затем настроить запись флагов в border-слой, чтобы
избежать отдельного использования функции Tourism::prepare.
God mode: обмениваться только короткими массивами изменившихся флагов.

Build_aliens надо оптимизировать с учетом того, что он строится на основе
существующих списков. Либо полностью изменить алгоритм пересылок внутри
функции apply_flags, тут надо думать.