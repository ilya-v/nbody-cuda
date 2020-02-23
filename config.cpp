#define _CRT_SECURE_NO_WARNINGS 1

#include "config.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>

#include "params.hpp"

double rand(long* rng_state);

void config_t::generate(const params_t& params, long& rng) {
    if (params.initial_config.empty()) {
        std::cout << "# Skipping generation of input.txt with initial config" 
                  << std::endl;
        return;
    }

    std::unique_ptr<FILE, int(*)(FILE*)> fout{ fopen("input.txt", "w"), 
                                               fclose };
    if (!fout) {
        std::cout << "# Error creating input.txt" << std::endl;
        return;
    }

    if (params.initial_config == "sphere") {
        std::cout << "# Generating initial configuration of type 'sphere'"
                  << std::endl;
        unsigned count = 0;
        const int Ns = (unsigned)(params.config_r / params.config_a) + 1;
        for (int i = -Ns; i <= Ns; i++)
        for (int j = -Ns; j <= Ns; j++)
        for (int k = -Ns; k <= Ns; k++) {
            const double
                x = params.config_a * i,
                y = params.config_a * j,
                z = params.config_a * k;
            if (x * x + y * y + z * z > params.config_r* params.config_r)
                continue;
            const double
                vx = -params.config_ang_vel * y,
                vy = params.config_ang_vel * x,
                vz = 0;
            const double
                vrx = 2 * (rand(&rng) - 0.5) * params.config_rand_vel,
                vry = 2 * (rand(&rng) - 0.5) * params.config_rand_vel,
                vrz = 2 * (rand(&rng) - 0.5) * params.config_rand_vel;

            fprintf(fout.get(), "%lf %lf %lf %lf %lf %lf\n",
                x, y, z, vx + vrx, vy + vry, vz + vrz);
            count++;
        }
        std::cout << "# Generated " << count << " particles" << std::endl;
    }
    else {
        std::cout << "# Unknown initial_config: <" << params.initial_config
            << ">" << std::endl;
    }
}


bool config_t::read(const std::string &fname,
                    std::function<void(std::array<double, 6>)> add_particle) {
    std::unique_ptr<FILE, int(*)(FILE*)> fin{ fopen(fname.c_str(), "r"),
                                              fclose };
    if (!fin) {
        std::cout <<  "# Error: cannot open input.txt" << std::endl;
        return false;
    }

    unsigned i = 0;
    int n_read = 0;
    for (;i == 0 || n_read == 6; i++) {
        std::array<double, 6> arr{};
        n_read = fscanf(fin.get(), "%lf  %lf  %lf  %lf  %lf  %lf\n",
                &arr[0], &arr[1], &arr[2], &arr[3], &arr[4], &arr[5]);
        if (6 == n_read)        
            add_particle(std::move(arr));
    }

    if (n_read != 0) {
        std::cout << "# Error reading line " << i << " from " << fname
            << ": only " << n_read << " values read." << std::endl;
        return false;
    }

    if (1 >= fscanf(fin.get(), "%*s")) {
        std::cout << "# Error: dangling data in " << fname <<" after reading "
                  << i << " lines." << std::endl;
        return false;
    }
    return true;
}
