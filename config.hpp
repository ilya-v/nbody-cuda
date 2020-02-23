#ifndef CONFIG__HPP_
#define CONFIG__HPP_

#include <array>
#include <functional>
#include <string>

struct params_t;

struct config_t {
    static void generate(const params_t& params, long& rng);
    static bool read(const std::string &fname, 
                     std::function<void(std::array<double, 6>)>);
};


#endif
