/// @brief В проекте реализована поддержка JSON, большинство классов
/// могут быть инициализированы из JSON строки.

#include <iostream>

#include <zephyr/utils/json.h>

using zephyr::utils::Json;

int main(int argc, char** argv) {
    Json config = Json::load("config.json");

    return 0;
}
