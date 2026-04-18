# Zephyr

> Современный высокопроизводительный CFD-фреймворк с локально-адаптивными сетками
> Документация доступна по ссылке: [https://saintluka.github.io/zephyr/](https://saintluka.github.io/zephyr/)

![CMake](https://img.shields.io/badge/build-CMake-brightgreen.svg)
![C++20](https://img.shields.io/badge/C++-20/23-blue.svg)
![MPI](https://img.shields.io/badge/MPI-enabled-important)
![TBB](https://img.shields.io/badge/TBB-parallel-important)

# Ключевые особенности

- Современный C++17/20: Чистый производительный код с использованием шаблонов и STL.
- Локально-адаптивные сетки: динамическое разбиение/укрупнение сетки (AMR) для повышения точности в критических областях.
- Высокопроизводительные вычисления: гибкий параллелизм с использованием многопоточности на узле и MPI для распределенных вычислений на кластерах.
- Кросс-платформенная сборка: Упрощенный процесс сборки с помощью CMake.
- Поддержка различных физических моделей:
  - Классическая газовая динамика
  - Многоматериальная динамика
- Модульная архитектура: Библиотека ядра для работы с сетками может использоваться независимо от решателей.

# Структура проекта

```
zephyr-src/
└── README.md          # Этот файл
├── CMakeLists.txt     # Главный файл конфигурации сборки
├── docs/              # Документация
├── problems/          # 
│   ├── single-fluid/  # Одноматериальные задачи
│   ├── multi-fluid/   # Многоматериальные задачи
├── tests/             # Модульные тесты
├── zephyr/            # Исходные файлы библиотеки
│   ├── geom/          # Модуль геометрии
│   ├── math/          # Математический модуль (решатели)
│   ├── mesh/          # Работа с сетками
│   └── utils/         # Многопоточность, поддержка MPI
```

# Предварительные требования

- Компилятор C++ с поддержкой C++20 (GCC >= 11, Clang >= 10, MSVC >= 2019)
- CMake >= 3.16
- MPI (опционально, для распределенных расчетов): OpenMPI, MPICH, Intel MPI
- Intel TBB (рекомендуется) или другой менеджер потоков

# Сборка и установка

Клонирование репозитория в `zephyr-src`:
```bash
git clone https://github.com/SaintLuka/zephyr.git zephyr-src
```

Установка зависимостей:
```bash
sudo apt install libeigen3-dev libboost-dev
```

Запомним корневую директорию:
```bash
export ROOT=$PWD
```

Создать директорию сборки и перейти в неё:
```bash
mkdir build && cd build
```

Конфигурация проекта:
```bash
cmake ../zephyr-src/ -DZEPHYR_PROBLEMS_ALL=ON -DCMAKE_INSTALL_PREFIX=$ROOT/install
```
*Совет*: Используйте ccmake или cmake-gui для интерактивной настройки.

Собрать проект и установить:
```bash
make -j8 install
```

Установить путь до библиотеки `libzephyr.so`:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOT/install/lib
```

Перейти в директорию с исполняемым файлом и запустить:
```bash
cd $ROOT/install/problems/single-fluid
./smf_test_2D --threads=on
```

Результаты выполнения программы появятся в папке output в формате VTK.
