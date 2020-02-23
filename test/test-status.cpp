#include "../status.hpp"

#include <cassert>
#include <iostream>

void test_status() {
    {
        status_t s = {0,};
        std::string h = s.header();
        std::string l = s.line();
        assert(l.size() == h.size());

    }
}