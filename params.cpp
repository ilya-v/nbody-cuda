#define _CRT_SECURE_NO_WARNINGS 1

#include "params.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include <cstdio>

const params_t params_t::default_params() {
    params_t p = { 0, };
    p.n_steps = 1000 * 1000;
    p.n_steps_for_report = 100;
    p.n_steps_for_output = 1000;
    p.dt_max = 0.01;
    p.dv_max = 0.01;
    p.dt_start = 0.001;
    p.t_start = 0;
    p.r2_eps = 0.01;
    p.initial_config = "sphere";
    p.config_a = 1.0;
    p.config_r = 10.0;
    p.config_ang_vel = 0.1;
    p.config_rand_vel = 0.1;
    p.random_seed = 102030;
    return p;
}

namespace {
    struct param_adapter_t {
        std::string name;
        std::function<bool(const std::string & line)> read;
        std::function<std::string()> print;

        template<typename T>
        param_adapter_t(const std::string& name, T& value) 
            : name{ name }, print{ printer(name, value) } {
            read = [name, &value](const std::string& line) mutable {
                char key[256] = { 0, }, strvalue[256] = { 0, };
                const char* li = line.c_str();
                if (2 == sscanf(li, " %[a-z_A-Z0-9] %*[=] %s", key, strvalue) &&
                    name == key) {
                    std::istringstream ss{strvalue};
                    ss >> value;
                    return !ss.fail();
                }
                return false;
            };
        }

        template<typename T>
        param_adapter_t(const std::string& name, const T& value)
            : name{ name }, print{printer(name, value)} {}

    private:
        template<typename T>
        std::function<std::string()> printer(const std::string& name, 
                                             const T& value) const {
            return [&value] {
                std::ostringstream ss;
                ss << value;
                return ss.str();
            };
        }
    };

    template<typename P>
    std::vector<param_adapter_t> param_adapters(P &&p) {
        return std::vector<param_adapter_t> {
            { "n_steps",            p.n_steps },
            { "n_steps_for_report", p.n_steps_for_report },
            { "n_steps_for_output", p.n_steps_for_output },
            { "dt_max",             p.dt_max },
            { "dv_max",             p.dv_max },
            { "dt_start",           p.dt_start },
            { "t_start",            p.t_start },
            { "r2_eps",             p.r2_eps },
            { "initial_config",     p.initial_config },
            { "config_a",           p.config_a },
            { "config_r",           p.config_r },
            { "config_ang_vel",     p.config_ang_vel },
            { "config_rand_vel",    p.config_rand_vel },
            { "random_seed",        p.random_seed}
        };
    }
}

void params_t::read(const std::string &fname) {
    std::ifstream fparam{fname};
    if (!fparam) {
        std::cout << "# " << fname 
                << " not found, using default parameters" 
                << std::endl;
        *this = default_params();
        return;
    }
    const std::vector<param_adapter_t> pas = param_adapters(*this);
    for (std::string line; std::getline(fparam, line);)
        for (const param_adapter_t& pa: pas)
            if (pa.read(line))
                break;
}

std::string params_t::to_string(const std::string &pre, 
                                const std::string &post) const {
    const std::vector<param_adapter_t> pas = param_adapters(*this);
    std::ostringstream ss;
    for (const param_adapter_t &pa : pas)
        ss  << std::left << pre << std::setw(20) 
            << pa.name << " = " << pa.print() << post;
    return ss.str();
}

