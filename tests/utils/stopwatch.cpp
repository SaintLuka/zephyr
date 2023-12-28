#include <iostream>
#include <thread>
#include <cmath>

#include <zephyr/utils/stopwatch.h>

using zephyr::utils::Stopwatch;

void sleep(int ms) {
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

int foo(int seed) {
    int res = seed;
    for (int i = 0; i < 30000000; ++i) {
        res += int(154214.0 * std::sin(res) * std::sin(res));
        res %= 113;
    }
    return res;
}

int main() {
    Stopwatch sw;

    sleep(100);

    std::cout << sw.milliseconds()<< " ms\n";

    sw.resume();
    sleep(40);
    sw.stop();

    std::cout << sw.milliseconds() << " ms\n";
    sleep(10);

    sw.resume();
    sleep(60);
    std::cout << sw.milliseconds() << " ms\n";
    sw.stop();

    sw.start();
    sleep(50);
    std::cout << sw.milliseconds() << " ms\n";
    sw.stop();

    sw.start();
    int a = foo(0);
    sw.stop();

    std::cout << "res " << a << "; " << sw.milliseconds() << " ms\n";

    int b = sw.measure(foo, 0);

    std::cout << "res " << b << "; " << sw.milliseconds() << " ms\n";

    return 0;
}