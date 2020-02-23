#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "config.hpp"
#include "params.hpp"
#include "status.hpp"

struct timer_t {
    clock_t t0 = clock();
    double time_sec() const {
        return (clock() - t0) / (double)CLOCKS_PER_SEC;
    }
};

#define BLOCK_SIZE 256

struct particle {
    double x, y, z, vx, vy, vz;
    particle() = default;
    particle(std::array<double, 6> arr) :
        x{ arr[0] }, y{ arr[1] }, z{ arr[2] },
        vx{ arr[3] }, vy{ arr[4] }, vz{ arr[5] }  {}
};

void status_update( status_t &status,
                    const unsigned N,
                    const particle* particles,
                    const std::vector<double> &u,
                    const double machine_time)
{
    double px = 0, py = 0, pz = 0;
    status.ek = 0;
    status.I = 0;
    for (const auto* p = particles; p < particles + N; p++) {
        px += p->vx; py += p->vy; pz += p->vz;
        status.ek += (p->vx*p->vx + p->vy*p->vy + p->vz*p->vz)/2;
        status.I += p->x*p->x + p->y*p->y + p->z*p->z;
    }

    status.ep = 0;
    for (double ui: u)
        status.ep += ui;

    status.etot = status.ek - status.ep;
    status.machine_time = machine_time;
}

void particles_center(std::vector<particle> &particles) {
    double  x = 0,  y = 0,  z = 0;
    double  px = 0, py = 0, pz = 0;
    for (const auto &p : particles) {
        px += p.vx; py += p.vy; pz += p.vz;
        x += p.x;   y += p.y;   z += p.z;
    }
    const auto N = particles.size();
    px/=N;  py/=N;  pz/=N;
    x/=N;   y/=N;   z/=N;
    for (auto &p : particles) {
        p.vx -= px; p.vy -= py; p.vz -= pz;
        p.x -= x;   p.y -= y;   p.z -= z;
    }
}

__global__
void calcForces(particle *p, double dt, unsigned N, double r2_eps) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double fx = 0, fy = 0, fz = 0;
        for (unsigned j = 0; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx * dx + dy * dy + dz * dz + r2_eps,
                invDist = rsqrt(distSqr),
                invDist3 = invDist * invDist * invDist;

            fx += dx * invDist3;
            fy += dy * invDist3;
            fz += dz * invDist3;
        }

        p[i].vx += dt*fx;
        p[i].vy += dt*fy;
        p[i].vz += dt*fz;
    }
}

__global__
void calcPotential(particle *p, double *u, unsigned N, double r2_eps) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double ui = 0;
        for (unsigned j = i + 1; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx * dx + dy * dy + dz * dz + r2_eps;
            ui += rsqrt(distSqr);
        }
        u[i] = ui;
    }
}

int main(const int argc, const char** argv) {

    params_t params = params_t::default_params();

    if (argc >= 2 && strstr(argv[1], "-p")) {
        std::cout << params.to_string("#", "\n");
        return 0;
    }
    params.read("params.txt");
    std::cout << params.to_string("#", "\n");

    status_t status{};
    std::cout << status.header() << std::endl;

    status.rng = -(long)params.random_seed;
    config_t::generate(params, status.rng);;

    std::vector<particle> particles;
    config_t::read("input.txt", [&particles](std::array<double, 6> arr) {
        particles.emplace_back(particle{std::move(arr)});
    });

    std::vector<particle> old_particles(particles.size(), particle{});
    assert(!particles.empty());
    particles_center(particles);

    const int N = static_cast<int>(particles.size());
    const int nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    particle *d_p;
    cudaMalloc(&d_p, N*sizeof(particle));

    double *d_u;
    cudaMalloc(&d_u, N*sizeof(double));

    status.t = params.t_start;
    status.dt = params.dt_start;

    particle *p_particles = &particles.front();
    particle *p_old_particles = &old_particles.front();

    timer_t timer{};
    for(status.step = 0; status.step <= params.n_steps; status.step++) {

        cudaMemcpy( d_p, p_particles, N*sizeof(particle),
                    cudaMemcpyHostToDevice);
        calcForces<<<nBlocks, BLOCK_SIZE>>>(d_p, status.dt, N, params.r2_eps);
        std::swap(p_particles, p_old_particles);
        cudaMemcpy( p_particles, d_p, N*sizeof(particle),
                    cudaMemcpyDeviceToHost);

        double dv = 0;
        for (int i = 0 ; i < N; i++) {
            particle *p = p_particles + i;
            p->x += p->vx*status.dt;
            p->y += p->vy*status.dt;
            p->z += p->vz*status.dt;

            particle *op = p_old_particles + i;
            const double
                dv_ix = p->vx - op->vx,
                dv_iy = p->vy - op->vy,
                dv_iz = p->vz - op->vz,
                dv_i = sqrt(dv_ix * dv_ix + dv_iy*dv_iy + dv_iz*dv_iz);
            if (dv_i > dv)
                dv = dv_i;
        }
        status.t += status.dt;

        status.dt = params.dv_max/dv * status.dt;
        if (status.dt > params.dt_max)
            status.dt = params.dt_max;

        if (status.step % params.n_steps_for_report == 0) {
            std::vector<double> u( static_cast<size_t>(N), 0.0 );
            cudaMemcpy( d_u, &u.front(), N*sizeof(double),
                        cudaMemcpyHostToDevice);
            calcPotential<<<nBlocks, BLOCK_SIZE>>>(d_p, d_u, N, params.r2_eps);
            cudaMemcpy(&u.front(), d_u, N*sizeof(double),
                        cudaMemcpyDeviceToHost);

            status_update(status, N, p_particles, u, timer.time_sec());
            std::cout << status.line() << std::endl;
        }

        if (status.step % params.n_steps_for_output == 0) {
            char fname[256];
            sprintf(fname, "out-%06u.xyz", status.step);
            FILE *fout = fopen(fname, "w");
            fprintf(fout, "%u\nt=%.8lf\n", N, status.t);
            for (int i = 0; i < N; i++) {
                particle *p = p_particles + i;
                fprintf(fout, "C %lf %lf %lf %lf %lf %lf\n",
                        p->x, p->y, p->z, p->vx, p->vy, p->vz);
            }
            fclose(fout);
        }
    }

    printf("#N=%d, Steps=%u Timer=%0.3lf s\n",
            N, params.n_steps, timer.time_sec() / params.n_steps);
    cudaFree(d_p);
    cudaFree(d_u);

    cudaDeviceSynchronize();
    cudaDeviceReset();
}
