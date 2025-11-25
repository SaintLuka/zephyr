---
sidebar_position: 4
---

# CMake

## Что такое CMake?


### Опции CMake


### CMakeCache


## Список опций

Полный список опций.

Включить/выключить компиляцию дополнительных библиотек можно с помощью опций CMake:
<div class="no-header-table">
| | | |
|-|-|-|
|`THREADS_TYPE` | STD / TBB | Использовать многопоточность из стандартной библиотеки (STD) или библиотеку TBB. |
|`ZEPHYR_MPI`    | ON / OFF | Включить компиляцию функций для распределенных вычислений. |
|`ZEPHYR_PYTHON` | ON / OFF | Включить библиотеки python, станут доступны средства визуализации. |
</div>

Следующие опции CMake контроллируют компиляцию исполняемых файлов:
<div class="no-header-table">
| | | |
|-|-|-|
|`ZEPHYR_TESTS`        | ON / OFF | Компиляция части тестов |
|`ZEPHYR_TESTS_ALL`    | ON / OFF | Компиляция всех тестов  |
|`ZEPHYR_PROBLEMS`     | ON / OFF | Компиляция части задач  |
|`ZEPHYR_PROBLEMS_ALL` | ON / OFF | Компиляция всех задач   |
</div>


Дополнительные опции, которые трогать не следует:
<div class="no-header-table">
| | | |
|-|-|-|
|`ZEPHYR_ASSERTS` | ON / OFF |
|`ZEPHYR_DOXYGEN` | ON / OFF |
|`ZEPHYR_EIGEN`   | ON / OFF |
</div>

