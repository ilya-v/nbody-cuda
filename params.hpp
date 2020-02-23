#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <vector>
#include <string>

struct params_t {
    unsigned
        n_steps,
        n_steps_for_report,
        n_steps_for_output;

    double
        dt_max,
        dv_max,
        dt_start,
        t_start,
        r2_eps;

    std::string initial_config;
    double
        config_a,
        config_r,
        config_ang_vel,
        config_rand_vel;

    unsigned random_seed;

    static const params_t default_params();

    void read(const std::string &fname);
    std::string to_string(const std::string &pre = "#",
                          const std::string &post = "\n") const;
};

#endif

