#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

#define BLOCK_SIZE 256
#define SOFTENING 1e-3f

typedef struct { double x, y, z, vx, vy, vz; } Particle;

__global__
void calcForces(Particle *p, double dt, int N) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

        for (int j = 0; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx*dx + dy*dy + dz*dz + SOFTENING,
                invDist = rsqrtf(distSqr),
                invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        p[i].vx += dt*Fx;
        p[i].vy += dt*Fy;
        p[i].vz += dt*Fz;
    }
    __syncthreads();
}

int main(const int argc, const char** argv) {

    const int
        Ns =  50,
        N  = Ns * Ns * Ns,
        nSteps = 1000,
        nStepsForReport = 10,
        nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    const double
        dt = 0.001f,
        init_dist = 1.0;

    size_t n_bytes = N*sizeof(Particle);
    Particle *particles = (Particle*)malloc(n_bytes);

    for (int i = 0; i < Ns; i++)
    for (int j = 0; j < Ns; j++)
    for (int k = 0; k < Ns; k++)  {
        const double init_shift = init_dist * (Ns - 1) / 2;
        Particle *p = particles + i;
        p->x = i * init_dist - init_shift;
        p->y = j * init_dist - init_shift;
        p->z = k * init_dist - init_shift;
        p->vx = 0;
        p->vy = 0;
        p->vz = 0;
    }

    Particle *d_p;
    cudaMalloc(&d_p, n_bytes);

    StartTimer();
    for (int iter = 0; iter < nSteps; iter++) {

        cudaMemcpy(d_p, particles, n_bytes, cudaMemcpyHostToDevice);
        calcForces <<<nBlocks, BLOCK_SIZE>>>(d_p, dt, N);
        cudaMemcpy(particles, d_p, n_bytes, cudaMemcpyDeviceToHost);

        for (int i = 0 ; i < N; i++) {
            particles[i].x += particles[i].vx*dt;
            particles[i].y += particles[i].vy*dt;
            particles[i].z += particles[i].vz*dt;
        }

        if (iter % nStepsForReport == 0) {
            double px = 0, py = 0, pz = 0;
            for (int i = 0 ; i < N; i++) {
                Particle *p = particles + i;
                px += p->vx;
                py += p->vy;
                pz += p->vz;
            }
            printf("p %lf %lf %lf\n", px, py, pz);
        }
    }

    printf("N=%d, Titer=%0.3lf s\n",
        N,
        GetTimer() / nSteps / 1000.0);

    free(particles);
    cudaFree(d_p);
}
