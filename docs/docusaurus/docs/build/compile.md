---
sidebar_position: 2
description: "Структура проекта, сборка на различных платформах"
toc_max_heading_level: 3  # Только h2 (первый уровень)
toc_min_heading_level: 2  # Начинать с h2
---

# Сборка и запуск

## Структура проекта

***Создание рабочей директории***

Для инициализации проекта необходимо создать корневую директорию и выполнить клонирование репозитория:

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

- **zephyr-src** – корневая директория исходного кода, точка подключения git;

- **build** – директория для сборки проекта;

- **zephyr** – исходные файлы библиотеки libzephyr;

- **problems** – исходные коды расчетных задач;

- **tests** – исходные коды тестов к модулям;

- **CMakeLists.txt** – основной конфигурационный файл системы сборки cmake.

Рекомендуется выполнять сборку в отдельной директории. При сборке внутри `zephyr-src` следует убедиться, что файлы компиляции не добавляются в git. Файлы компиляции, а также файлы конфигураций, созданные IDE, следует добавить в файл .gitignore.

## Поддерживаемые платформы

1. **Windows (MSVC)**. Базовая конфигурация с ограничениямя, вариант подходит для быстрого старта. Есть сложности с проведением распределенных вычислений (MPI) и интеграцией библиотек на python.

2. **Windows (MinGW)**. ?

3. **Linux (GCC/CLang)**. Рекомендуемая платформа, полная функциональность. Самый простой вариант настройки зависимостей, поскольку систему сборки CMake можно считать родной для мира open source.

## Windows (MSVC)

<details style={{
  border: "none", boxShadow: "none",
  background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### Требования к окружению

- Операционная система Windows 10/11;
- Компилятор Visual Studio Build Tools 2022/2026;
- Поддержка стандарта C++20.

Есть два варианта установки: полностью установить Visual Studio, или только Visual Studio Build Tools (компилятор), если по какой-то причине нет желания использовать Visual Studio. В любом случае запускаем установщик и выбираем пресет "Разработка классических приложений на C++". Отмечаем пункты:
- C++ Build Tools
- Windows SDK
- Средства CMake C++

Рекомендую не включать пункт "Диспетчер пакетов vcpkg", а установить vcpkg самостоятельно. После установки Build Tools станет доступна командная строка "x64 Native Tools Command Prompt for VS", далее все команды выполняются в ней.

### Установка зависимостей

Для работы с зависимостями используем диспетчер пакетов **vcpkg**. Клонируем репозиторий в корень C:\\
```bash
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
bootstrap-vcpkg.bat
```

В папке **C:\vcpkg** есть исполняемый файл vcpkg, который используется для установки обязательных библиотек:
```bash
vcpkg install boost
vcpkg install eigen3
vcpkg install tbb
vcpkg integrate install
```
Boost устанавливается довольно долго (до двух часов). Вызов `integrate install` выдает путь к конфигурации CMake, который следует запомнить. В данном случае это

`-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake`

Библиотека TBB необязательная, но у меня не возникло проблем с установкой.

### Сборка и запуск

***Конфигурация CMake***

Перейдем в консоли в директорию сборки. Выполним конфигурацию CMake:
```bash
mkdir build && cd build
cmake ..\zephyr-src -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
```

Дополнительные ключи позволяют конкретизировать компилятор и платформу:
```bash
cmake ..\zephyr-src -G "Visual Studio 17 2022" -A x64
```

При работе из "x64 Native Tools Command Prompt for VS" в этом нет необходимости. Дополнительные опции CMake (build type, install prefix и другие) можно передать при первом вызове, а можно запустить конфигурацию CMake несколько раз, добавляя или изменяя опции. К примеру, включим сборку всех задач и используем библиотеку TBB:
```bash
cmake ..\zephyr-src -DTHREADS_TYPE=TBB
cmake ..\zephyr-src -DZEPHYR_PROBLEMS_ALL=ON
```

<details>
<summary>Конфигурация с помощью cmake-gui</summary>

После запуска программы:
1. Выбрать папку с исходниками (с основным CMakeLists.txt);
2. Выбрать папку для сборки (build);
3. Выставить toolset;
4. Очистить кэш (File -> Delete Cache);
5. Добавить переменную (Add entry) `CMAKE_TOOLCHAIN_FILE: "C:\vcpkg\scripts\buildsystems\vcpkg.cmake"`

Теперь можно визуально увидеть зависимости, добавить или изменить опции CMake, вручную прописать пути к библиотекам и так далее.

</details>


***Сборка проекта***

Сборка в режиме Debug:

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

***Запуск в коносли***

Библиотека libzephyr.dll является динамической, поэтому при запуске исполняемый файл должен "знать", где она находится.
Библиотеку можно скопировать и положить рядом с исполняемым файлом, также можно создать символическую ссылку на библиотеку рядом с исполняемым файлом. А можно добавить путь к библиотке в переменную окружения `PATH`. Выберем последний вариант и запустим задачу `wave`:

```bash
set PATH=C:\...path...\build\bin\Release;%PATH%
cd problmes\Release
wave
```

Результаты выполнения появятся рядом с программой в папке output. Для визуализации следует использовать
[ParaView](/docs/additional/paraview).


<details style={{
  border: "none", boxShadow: "none",
  background: "none", padding: "0"}}>
<summary>
**Нерешенные проблемы**
</summary>

Проблемы с установкой на некоторых платформах, которые пока не решены. Будет классно, если кто-то их починит.

***Использование встроенного vcpkg***

Менеджер пакетов vcpkg, встроенный в Visual Studio работает иначе. Для него невозможна глобальная установка пакетов. Делается по аналогии с виртуальным окружением для python, то есть зависимости устанавливаются локально рядом с проектом.

Если вызвать vcpkg integrate install,
то также выплюнет путь.

В папке проекта выполняем

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


## Windows (MinGW)

<details style={{
  border: "none", boxShadow: "none",
  background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### Компиляторы и средства разработки

### Установка зависимостей

pacman?

### Сборка и запуск

</details>


## Linux (GCC/Clang)

<details style={{
  border: "none", boxShadow: "none",
  background: "none", padding: "0"}}>
<summary>Инструкция</summary>

### Установка зависимостей

### Сборка и запуск

</details>


