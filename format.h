#pragma once
#if __has_include(<format>)
#include <format>
#else
#include <string>
namespace std {
  template <typename ... Args>
  std::string format(Args... args) {
    return "has no std::format, sorry";
  }
}
#endif