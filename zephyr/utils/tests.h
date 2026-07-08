#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <iomanip>
#include <vector>
#include <functional>
#include <stdexcept>
#include <format.h>

// ========== ИСКЛЮЧЕНИЕ ДЛЯ ТЕСТОВ ==========
class TestFailure : public std::runtime_error {
public:
    explicit TestFailure(const std::string& msg) : std::runtime_error(msg) {}
};

// ============ ОСНОВНЫЕ МАКРОСЫ =============

/// @brief Проверка условия
#define TEST_ASSERT(expr) \
    do { \
        if (!(expr)) { \
            throw TestFailure( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": Assertion failed: " #expr \
            ); \
        } \
    } while(0)

/// @brief Проверка на равенство
#define TEST_ASSERT_EQ(a, b) \
    do { \
        auto _a = (a); \
        auto _b = (b); \
        if (!(_a == _b)) { \
            throw TestFailure( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": " #a " == " #b " failed (" + std::to_string(_a) + " != " + std::to_string(_b) + ")" \
            ); \
        } \
    } while(0)

/// @brief Совпадает в шести знаках
#define TEST_ASSERT_CLOSE(a, b) \
    do { \
        auto _a = (a); \
        auto _b = (b); \
        double _eps = 1.0e-6 * std::max(std::abs(_a), std::abs(_b)); \
        if (!(std::abs(_a - _b) <= _eps)) { \
            throw TestFailure(std::format("{}:{}: |{} - {}| <= eps failed (|{:.3e} - {:.3e}| > {:.3e})", \
                             __FILE__, __LINE__, #a, #b, _a, _b, _eps)); \
        } \
    } while(0)

/// @brief Проверка на неравенство
#define TEST_ASSERT_NE(a, b) \
    do { \
        auto _a = (a); \
        auto _b = (b); \
        if (!(_a != _b)) { \
            throw TestFailure( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": " #a " != " #b " failed (" + std::to_string(_a) + " == " + std::to_string(_b) + ")" \
            ); \
        } \
    } while(0)

/// @brief Проверка строк на равенство
#define TEST_ASSERT_STR_EQ(a, b) \
    do { \
        std::string _a = (a); \
        std::string _b = (b); \
        if (_a != _b) { \
            throw TestFailure( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": Strings not equal: \"" + _a + "\" != \"" + _b + "\"" \
            ); \
        } \
    } while(0)

/// @brief Проверка на вызов исключения
#define TEST_ASSERT_THROWS(expr) \
    do { \
        bool _caught = false; \
        try { expr; } \
        catch (...) { _caught = true; } \
        if (!_caught) { \
            throw TestFailure( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": Expected exception, but none thrown: " #expr \
            ); \
        } \
    } while(0)

// ========== РЕГИСТРАЦИЯ ТЕСТОВ ==========
class TestRegistry {
public:
    using TestFunc = std::function<void()>;

    static TestRegistry& instance() {
        static TestRegistry registry;
        return registry;
    }

    void add(const std::string& name, TestFunc func) {
        m_tests.push_back({name, std::move(func)});
    }

    int run() const {
        int passed = 0;
        int failed = 0;

        std::cout << "Running " << m_tests.size() << " test(s)...\n";

        for (const auto& [name, func] : m_tests) {
            std::cout << "  " << name << "... ";
            try {
                func();
                std::cout << "OK\n";
                passed++;
            } catch (const TestFailure& e) {
                std::cout << "FAIL\n";
                std::cerr << "    " << e.what() << "\n\n";
                failed++;
            } catch (const std::exception& e) {
                std::cout << "ERROR\n";
                std::cerr << "    Unexpected exception: " << e.what() << "\n\n";
                failed++;
            } catch (...) {
                std::cout << "ERROR\n";
                std::cerr << "    Unknown exception\n\n";
                failed++;
            }
        }

        std::cout << "\n----------- Results -----------\n";
        std::cout << "  Passed: " << passed << "\n";
        std::cout << "  Failed: " << failed << "\n";
        std::cout << "  Total:  " << passed + failed << "\n";
        std::cout << "-------------------------------\n";

        return failed == 0 ? 0 : 1;
    }

private:
    struct TestCase {
        std::string name;
        TestFunc func;
    };

    std::vector<TestCase> m_tests;
};

// =============== МАКРОС ДЛЯ ТЕСТОВ =============
#define CONCAT_IMPL(x, y) x##y
#define CONCAT(x, y) CONCAT_IMPL(x, y)

#define TEST_CASE_IMPL(name, id) \
    static void CONCAT(test_func_, id)(); \
    namespace { \
        struct CONCAT(TestRegistrar_, id) { \
            CONCAT(TestRegistrar_, id)() { \
                TestRegistry::instance().add(name, CONCAT(test_func_, id)); \
            } \
        }; \
        static CONCAT(TestRegistrar_, id) CONCAT(registrar_, id); \
    } \
    static void CONCAT(test_func_, id)()

#define TEST_CASE(name) TEST_CASE_IMPL(name, __LINE__)

