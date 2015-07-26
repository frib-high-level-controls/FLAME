#ifndef UTIL_H
#define UTIL_H

#include <stdexcept>

struct key_error : public std::runtime_error
{
    key_error(const std::string& s) : std::runtime_error(s) {}
};

#endif // UTIL_H
