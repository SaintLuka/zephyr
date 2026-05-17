#!/bin/bash

# Run script in build directory, like
# ./../zephyr-src/compile-check.sh ../zephyr-src

if [ $# -eq 0 ]; then
    echo "Usage: $0 <path_to_source_dir>"
    echo "Example: $0 /path/to/your/cmake_project"
    exit 1
fi

SOURCE_DIR="$(cd "$1" && pwd)"
CURRENT_DIR="$(pwd)"

# TODO: уменьшить список вариантов
CMAKE_OPTIONS=(
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=OFF -DZEPHYR_PROBLEMS_ALL=OFF"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=OFF -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=OFF -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=OFF -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
    "-DZEPHYR_ASSERTS=ON  -DZEPHYR_MPI=ON  -DZEPHYR_PYTHON=ON  -DZEPHYR_TESTS_ALL=ON -DZEPHYR_PROBLEMS_ALL=ON"
)

clean_options_string() {
    local input="$1"
    result=$(echo "$input" | sed 's/[^a-zA-Z0-9]/_/g')
    result="${result//DZEPHYR/}"
    result="${result//__/_}"
    echo "$result"
}

build_with_options() {
    local options="$1"
    local build_name=$(clean_options_string "$options")
    local build_dir="${CURRENT_DIR}/build_${build_name}"

    echo "========================================="
    echo "Testing with options: $options"
    echo "Build directory: $build_dir"
    echo "========================================="

    mkdir -p "$build_dir"
    cd "$build_dir" || return 1

    echo "Running CMake with source dir: $SOURCE_DIR"
    cmake "$SOURCE_DIR" $options || {
        echo "❌ CMake failed with options: $options"
        cd "$CURRENT_DIR"
        return 1
    }

    echo "Running make..."
    make -j$(nproc) || {
        echo "❌ Make failed with options: $options"
        cd "$CURRENT_DIR"
        return 1
    }

    echo "✅ Success with options: $options"
    cd "$CURRENT_DIR"
    return 0
}

main() {
    local failed=0

    for options in "${CMAKE_OPTIONS[@]}"; do
        build_with_options "$options"
        if [ $? -ne 0 ]; then
            failed=1
            break
        fi
    done

    if [ $failed -eq 0 ]; then
        echo "========================================="
        echo "✅ All builds completed successfully!"
    else
        echo "========================================="
        echo "❌ Build failed, stopping..."
        exit 1
    fi
}

main "$@"