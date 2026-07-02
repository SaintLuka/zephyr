/// @brief
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>

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
        tests_.push_back({name, func});
    }
    
    int run_all() const {
        int passed = 0;
        int failed = 0;

        std::cout << "\nRunning " << tests_.size() << " test(s)...\n\n";

        for (const auto& test : tests_) {
            std::cout << "  " << test.name << "... ";
            try {
                test.func();
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

        std::cout << "\n========== RESULTS ==========\n";
        std::cout << "Passed: " << passed << "\n";
        std::cout << "Failed: " << failed << "\n";
        std::cout << "Total:  " << passed + failed << "\n";
        std::cout << "=============================\n";

        return failed == 0 ? 0 : 1;
    }

private:
    struct TestCase {
        std::string name;
        TestFunc func;
    };

    std::vector<TestCase> tests_;
};

// ========== МАКРОС ДЛЯ ТЕСТОВ (С __COUNTER__) ==========
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

