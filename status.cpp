#define _CRT_SECURE_NO_WARNINGS 1
#include "status.hpp"

#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace {

    struct status_field_t {
        template <typename T>
        status_field_t(const T status_t::* p, const char* name, const char* fmt)
            : name{ name }, print_f{ [p, fmt](const status_t& s) {
                char buf[256] = {};
                sprintf(buf, fmt, s.*p);
                return std::string{ buf };
            } } {}

            std::string value(const status_t& s) const { return print_f(s); }

            std::string header() const {
                char buf[256] = {};
                sprintf(buf, "%*s", (int)(print_f(status_t{}).size()), name);
                return buf;
            }

    private:
        using print_func_t = std::function<std::string(const status_t&)>;
        const char* const name;
        const print_func_t print_f;
    };

    std::vector<status_field_t> status_fields{
        {   &status_t::step,        "step",     "%8u"       },
        {   &status_t::t,           "time",     "%16.8lf"   },
        {   &status_t::dt,          "dt",       "%16.8le"   },
        {   &status_t::ep,          "ep",       "%16.8le"   },
        {   &status_t::ek,          "ek",       "%16.8le"   },
        {   &status_t::etot,        "etot",     "%16.8le"   },
        {   &status_t::I,           "I",        "%16.8le"   },
        {   &status_t::machine_time,"mtime",    "%16.8lf"   },
        {   &status_t::rng,         "rng",      "%8u"       },
    };
}

std::string status_t::header() const {
    std::string result;
    for (const status_field_t& sf : status_fields)
        result += sf.header();
    return result;
}

std::string status_t::line() const {
    std::string result;
    for (const status_field_t& sf : status_fields)
        result += sf.value(*this);
    return result;
}