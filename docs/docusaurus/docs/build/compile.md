---
sidebar_position: 2
description: "Структура проекта, сборка на различных платформах"
toc_max_heading_level: 3  # Только h2 (первый уровень)
toc_min_heading_level: 2  # Начинать с h2
---

# Сборка и запуск

## Структура проекта

***Создание рабочей директории***

Создаем корневую директорию и выполняем клонирование репозитория:

```bash
mkdir Zephyr && cd Zephyr
git clone https://github.com/SaintLuka/zephyr.git zephyr-src
```

***Рекомендуемая структура проекта:***

```
Zephyr/
 ├─ build/
 └─ zephyr-src/
     ├─ .git/
     ├─ docs/
     ├─ problems/
     ├─ tests/
     ├─ zephyr/
     ├─ .gitignore
     ├─ CMakeLists.txt
     └─ README.md
```

***Описание компонентов:***

- **build** – директория для сборки проекта;

- **zephyr-src** – корневая директория исходного кода, точка подключения git;

- **zephyr** – исходные файлы библиотеки libzephyr;

- **problems** – исходные коды расчетных задач;

- **tests** – исходные коды тестов к модулям;

- **CMakeLists.txt** – основной конфигурационный файл системы сборки cmake.

Не рекомендуется выполнять сборку в директории с исходными файлами. Также следует убедиться, что файлы компиляции и файлы конфигураций, созданные IDE, не добавляются в git (для исключения файлов используется *.gitignore*).

## Поддерживаемые платформы

1. **Windows (MSVC)**. Сборка нативным компилятором Windows. Есть ограничения функциональности: проблемы с проведением распределенных вычислений (MPI) и интеграцией библиотек на python.

2. **Windows (UCRT)**. Компиляция в Windows с поддержкой UNIX библиотек и компиляторов. Есть ограничения функциональности: проблемы с проведением распределенных вычислений (MPI) и интеграцией библиотек на python.

3. **Linux (GCC/CLang)**. Рекомендуемая платформа, полная функциональность. Самый простой вариант настройки зависимостей, поскольку систему сборки CMake можно считать стандартом в разработке свободного программного обеспечения.

## Windows (MSVC)

<details style={{border: "none", boxShadow: "none", background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### Требования к окружению

- Операционная система Windows 10/11;
- Компилятор Visual Studio Build Tools 2022/2026;
- Поддержка стандарта C++20.

Варианты установки: полная установка Visual Studio, установка только компилятора Visual Studio Build Tools, если по какой-то причине нет желания использовать Visual Studio. В любом случае после запуска установщика выбираем пресет "Разработка классических приложений на C++". Отмечаем пункты:
- C++ Build Tools
- Windows SDK
- Средства CMake C++

Рекомендуется не включать пункт "Диспетчер пакетов vcpkg", а установить vcpkg самостоятельно. После установки Build Tools станет доступен терминал "x64 Native Tools Command Prompt for VS", далее все команды выполняются в нём.

### Установка зависимостей

Для работы с зависимостями используем диспетчер пакетов **vcpkg**. Клонируем репозиторий в корень C:\\
```bash
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
bootstrap-vcpkg.bat
```

В директории **C:\vcpkg** находится исполняемый файл vcpkg, который используется для установки обязательных библиотек:
```bash
vcpkg install boost
vcpkg install eigen3
vcpkg install tbb
vcpkg integrate install
```

Установка boost может длиться до двух часов. Библиотека TBB необязательна, но обычно установка не вызывает проблем. Вызов `integrate install` выводит путь к конфигурации CMake, который следует запомнить. Строка имеет следующий вид:

`-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake`

### Сборка и запуск

***Конфигурация CMake***

Переходим в директорию сборки и выполняем конфигурацию CMake:
```bash
mkdir build && cd build
cmake ..\zephyr-src -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
```

Дополнительные ключи позволяют указать компилятор и платформу:
```bash
cmake ..\zephyr-src -G "Visual Studio 17 2022" -A x64
```

При работе из терминала "x64 Native Tools Command Prompt for VS" в этом нет необходимости. Дополнительные опции CMake (build type, install prefix и другие) можно передать при первом вызове, а можно запустить конфигурацию CMake несколько раз, добавляя или изменяя опции. К примеру, включим сборку всех задач и используем библиотеку TBB:
```bash
cmake ..\zephyr-src -DTHREADS_TYPE=TBB
cmake ..\zephyr-src -DZEPHYR_PROBLEMS_ALL=ON
```

<details>
<summary>Конфигурация с помощью cmake-gui</summary>

Программа cmake-gui позволяет интерактивно настраивать опции CMake. После запуска cmake-gui необходимо:
1. Выбрать папку с исходниками (с основным CMakeLists.txt);
2. Выбрать папку для сборки (build);
3. Выставить toolset;
4. Очистить кэш (File -> Delete Cache);
5. Добавить переменную (Add entry) `CMAKE_TOOLCHAIN_FILE: "C:\vcpkg\scripts\buildsystems\vcpkg.cmake"`

Теперь можно визуально увидеть зависимости, добавить или изменить опции CMake, вручную прописать пути к библиотекам и так далее.

</details>


***Сборка проекта***

Сборка в режиме Debug (*по умолчанию*):

```bash
cmake --build . --config Debug
```

Сборка в режиме Release в четырёх потоках:

```bash
cmake --build . --config Release -j 4
```

Сборка конкретной цели:

```bash
cmake --build . --target wave --config Release
```

Если не указать Release, то сборка по умолчанию в Debug.

***Запуск в терминале***

Библиотека libzephyr.dll является динамической, поэтому при запуске исполняемый файл должен "знать", где она находится.
Библиотеку можно скопировать и положить рядом с исполняемым файлом, также можно создать символическую ссылку на библиотеку рядом с исполняемым файлом. А можно добавить путь к библиотеке в переменную окружения `PATH`. Выберем последний вариант и запустим задачу `wave`:

```bash
set PATH=C:\...path...\build\bin\Release;%PATH%
cd problmes\Release
wave
```

Результаты выполнения появятся рядом с программой в папке output. Для визуализации следует использовать
[ParaView](/docs/additional/paraview).


<details style={{border: "none", boxShadow: "none", background: "none", padding: "0"}}>
<summary>
**Нерешенные проблемы**
</summary>

Проблемы с установкой на некоторых платформах, которые пока не решены. Будет классно, если кто-то их починит.

***Использование встроенного vcpkg***

Менеджер пакетов vcpkg, встроенный в Visual Studio работает иначе. Для него невозможна глобальная установка пакетов. Делается по аналогии с виртуальным окружением для python, то есть зависимости устанавливаются локально рядом с проектом.

Если вызвать vcpkg integrate install,
то также выплюнет путь.

В директории проекта выполняем

```bash
vcpkg new --application
```

После этого перечислить зависимости в vcpkg.json

```
{
  "dependencies": ["boost", "eigen3", ...]
}
```

и выполнить vcpkg install. Здесь появились какие-то проблемы, `vcpkg install` у меня не сработало.

***Использование python в vcpkg***

Устанавливает легко
Дополнительные зависимости
```bash
vcpkg install pybind11
```

Внутрь vcpkg также установится python, на него можно установить дополнительные пакеты через pip.
Домашняя директория `PYTHONHOME=C:\vcpkg\installed\x64-windows\tools\python3`

Переходим и устанавливаем
```bash
cd C:\vcpkg\installed\x64-windows\tools\python3
python -m ensurepip --upgrade
python -m pip install --upgrade pip
```

Создаст папку Scripts, там можно вызвать pip/pip3
```bash
cd Scripts
pip3 install numpy
pip3 install matplotlib
```

Один раз удалось запустить почти интерактивный python, и то отображение через браузер. Без этого всё не нужно.

***Параллельность MPI с MSVC***

Наивная установка MsMpi через vcpkg ни разу не сработала
```bash
vcpkg install msmpi
```
Другие варианты установки или MPICH я не пробовал.
</details>

</details>


## Windows (UCRT)

<details style={{border: "none", boxShadow: "none",background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### MSYS2

Рекомендуется использовать платформу **MSYS2**.
MSYS2 — это программная платформа для Windows, которая предоставляет среду, эмулирующую Unix-подобные системы (в частности, Linux).

Скачать MSYS2 можно с официального сайта https://www.msys2.org/.

После установки станет доступен терминал Msys2, а также несколько специализированных терминалов (окружений): UCRT64, MinGW64, Clang. Во всех терминалах используются Unix-команды. Рекомендуется использовать терминал **UCRT64**, все последующие команды выполняются в этом терминале.

### Установка зависимостей

MSYS предоставляет менеджер пакетов для С и С++ **pacman**.

После первого запуска терминала UCRT выполним обновление всех пакетов:
```bash
pacman -Syu
```
После **перезагрузки** терминала снова выполняем
```bash
pacman -Su
```

Устанавливаем средства разработки (make, cmake):
```bash
pacman -S mingw-w64-ucrt-x86_64-toolchain
pacman -S mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-make
```

Установка зависимостей:
```bash
pacman -S mingw-w64-ucrt-x86_64-eigen3
pacman -S mingw-w64-ucrt-x86_64-boost
pacman -S mingw-w64-ucrt-x86_64-boost-libs
pacman -S mingw-w64-ucrt-x86_64-tbb
```

Библиотека TBB необязательна, но обычно установка не вызывает проблем.

### Сборка и запуск

***Конфигурация CMake***

В терминалах используется Unix-подобный синтаксис. В путях используется прямой слеш, а корневая директория обозначается **/C/**. Перейдем в директорию сборки:
```bash
cd /C/Users/.../Zephyr/build
```

Конфигурация CMake:
```bash
cmake ../zephyr-src/ -DCMAKE_BUILD_TYPE=Release
```

Опции CMake (build type, install prefix и другие) можно передать при первом вызове, а можно запустить конфигурацию CMake несколько раз, добавляя или изменяя опции. К примеру, включим сборку всех задач и используем библиотеку TBB:
```bash
cmake ../zephyr-src -DTHREADS_TYPE=TBB
cmake ../zephyr-src -DZEPHYR_PROBLEMS_ALL=ON
```

***Сборка проекта***

Сборка в режиме Debug (*по умолчанию*):

```bash
cmake --build . --config Debug
```

Сборка в режиме Release в четырёх потоках:

```bash
cmake --build . --config Release -j 4
```

Сборка конкретной цели:

```bash
cmake --build . --target wave --config Release
```

Если не указать Release, то сборка по умолчанию в Debug.

***Запуск в терминале***

Библиотека libzephyr.dll является динамической, поэтому при запуске исполняемый файл должен "знать", где она находится.
Добавим путь к библиотеке в переменную окружения `PATH`:
```bash
export PATH="/C/Users/.../Zephyr/build:$PATH"
```
*Снова видим смесь Unix и Windows: используется Windows переменная окружения PATH, но она задается в Unix стиле. В частности, используется разделитель ':', в то время как в переменной окружения PATH пути разделяются символом ';'.*

Запустим задачу `wave`:
```bash
cd problmes
./wave
```

Результаты выполнения появятся рядом с программой в папке output. Для визуализации следует использовать
[ParaView](/docs/additional/paraview).

<details>
<summary>
Переменные окружения
</summary>
Программы, собранные с UCRT, можно запускать из обычного терминала или двойным нажатием кнопки мыши. Для этого необходимо добавить в переменную окружения PATH пути до библиотек ucrt64 и libzephyr.
Для этого необходимо:

- Найти на панели Windows "Изменение системных переменных среды".

- Установить "PATH=\C:\msys64\ucrt64\bin;C:\Users\...\Zephyr\build;".

Важно поставить ';' в конце! Данным способом переменная устанавливается навсегда, после этого исполняемые файлы будут всегда запускаться. Минус такого подхода в том, что путь будет всегда указывать на библиотеку в Release сборке, поэтому для Debug тестов придется либо менять переменную, либо делать Debug сборку в той же директории.
</details>

<details style={{border: "none", boxShadow: "none", background: "none", padding: "0"}}>
<summary>
**Нерешенные вопросы**
</summary>

Запустить MPI и python.

</details>

</details>


## Linux (GCC/Clang)

<details style={{
  border: "none", boxShadow: "none",
  background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### Требования к окружению

Для сборки проекта необходим компилятор **GCC** версии **11.0** или выше с поддержкой стандарта C++20.
В современных дистрибутивах (например, Ubuntu 22.04 LTS и свежее) необходимые инструменты доступны в репозиториях по умолчанию и устанавливаются через стандартный менеджер пакетов apt.

### Установка зависимостей

Установка необходимых пакетов:
```bash
sudo apt update
sudo apt install libeigen3-dev libboost-all-dev
```

Установка необязательных пакетов (см. [зависимости](/docs/build/dependencies#используемые-библиотеки)):
```bash
sudo apt install libtbb-dev libopenmpi-dev pybind11-dev
sudo apt install python3-numpy python3-matplotlib
```

### Сборка и запуск

***Конфигурация CMake***

Переходим в директорию сборки, выполняем конфигурацию CMake:
```bash
mkdir build && cd build
cmake ../zephyr-src -DCMAKE_BUILD_TYPE=Release
```

Дополнительные опции CMake (build type, install prefix и другие) можно передать при первом вызове, а можно запустить конфигурацию CMake несколько раз, добавляя или изменяя опции.

К примеру, включим сборку всех задач и используем библиотеку TBB:
```bash
cmake ../zephyr-src -DTHREADS_TYPE=TBB
cmake ../zephyr-src -DZEPHYR_PROBLEMS_ALL=ON
```

<details>
<summary>Конфигурация с помощью cmake-curses-gui</summary>

**cmake-curses-gui** это консольная версия CMake с псевдографическим интерфейсом для интерактивной настройки проекта в консоли. Установка:
```bash
sudo apt install cmake-curses-gui
```

Теперь вместо команды cmake вводим:
```bash
ccmake ../zephyr-src
```

Теперь в терминале можно выбрать опции, тип сборки, указать путь для установки.
После завершения настройки нажать 'c' (configure), затем 'g' (generate).
Использование cmake-curses-gui очень удобно, особенно при сборке на кластере.

</details>


***Сборка проекта***

Сборка всего проекта (*ключ -j указывает число потоков*):

```bash
make -j 4
```

Сборка конкретной цели:
```bash
cmake -j 4 wave
```

***Запуск в терминале***

Перейдем в папку с задачами и запустим `wave`:
```bash
cd problmes
./wave
```

Результаты выполнения появятся рядом с программой в папке output. Для визуализации следует использовать [ParaView](/docs/additional/paraview).

</details>
