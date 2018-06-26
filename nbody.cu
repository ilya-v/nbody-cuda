#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static clock_t t0;

void startTimer() {
    t0 = clock();
}

double getTimer() {
    clock_t t = clock();
    return (t - t0) / (double)CLOCKS_PER_SEC;
}

#define BLOCK_SIZE 256
#define SOFTENING (0.01)

typedef struct { double x, y, z, vx, vy, vz; } Particle;

__global__ void calcForces(Particle *p, double dt, unsigned N) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double fx = 0, fy = 0, fz = 0;
        for (unsigned j = 0; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx*dx + dy*dy + dz*dz + SOFTENING,
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

__device__ double d_potential = 0;

__global__ void calcPotential(Particle *p, double *u, unsigned N) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double ui = 0;
        for (unsigned j = i + 1; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
            ui += rsqrt(distSqr);
        }
        u[i] = ui;
    }
}


int main(const int argc, const char** argv) {

    const int
        nSteps = 10000,
        nStepsForReport = 10,
        nStepsForOutput = nSteps/10;

    const double
        dt_max = 0.01,
        dv_max = 0.01;
    double
        dt = dt_max/10;

    Particle
        *particles = NULL,
        *old_particles = NULL;
    double *u = NULL;
    unsigned N = 0;

    {
        FILE *fin = fopen("input.txt", "rb");
        if (!fin) {
            printf("Cannot open input.txt\n");
            return -1;
        }

        {
            double t1, t2, t3, t4, t5, t6;
            for (; 6 == fscanf(fin, "%lf  %lf  %lf  %lf  %lf  %lf\n",
                               &t1, &t2, &t3, &t4, &t5, &t6);
                   N++);
        }

        rewind(fin);
        particles = (Particle *) malloc(N*sizeof(Particle));
        old_particles = (Particle *) malloc(N*sizeof(Particle));
        unsigned i = 0;
        for (; i < N; i++)
            if (6 != fscanf(fin, "%lf  %lf  %lf  %lf  %lf  %lf\n",
                                    &particles[i].x,
                                    &particles[i].y,
                                    &particles[i].z,
                                    &particles[i].vx,
                                    &particles[i].vy,
                                    &particles[i].vz))
                break;

        fclose(fin);

        if (i < N) {
            printf("Cannot read input.txt: %u lines from %u\n", i, N);
            return -1;
        }
        u = (double *) malloc(N * sizeof(double));
    }

    {
        double px = 0, py = 0, pz = 0;
        for (unsigned i = 0; i < N; i++) {
            Particle *p = particles + i;
            px += p->vx;
            py += p->vy;
            pz += p->vz;
        }
        px/=N;
        py/=N;
        pz/=N;
        for (unsigned i = 0; i < N; i++) {
            Particle *p = particles + i;
            p->vx -= px;
            p->vy -= py;
            p->vz -= pz;
        }

    }

    const int nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    Particle *d_p;
    cudaMalloc(&d_p, N*sizeof(Particle));

    double *d_u;
    cudaMalloc(&d_u, N*sizeof(double));

    double t = 0;
    startTimer();
    for (unsigned step = 0; step < nSteps; step++) {

        cudaMemcpy(d_p, particles, N*sizeof(Particle), cudaMemcpyHostToDevice);
        calcForces<<<nBlocks, BLOCK_SIZE>>>(d_p, dt, N);
        Particle *tmp_particles = old_particles;
        old_particles = particles;
        particles = tmp_particles;
        cudaMemcpy(particles, d_p, N*sizeof(Particle), cudaMemcpyDeviceToHost);

        double dv = 0;
        for (int i = 0 ; i < N; i++) {
            Particle *p = particles + i;
            p->x += p->vx*dt;
            p->y += p->vy*dt;
            p->z += p->vz*dt;

            Particle *op = old_particles + i;
            const double
                dv_ix = p->vx - op->vx,
                dv_iy = p->vy - op->vy,
                dv_iz = p->vz - op->vz,
                dv_i = sqrt(dv_ix * dv_ix + dv_iy*dv_iy + dv_iz*dv_iz);
            if (dv_i > dv)
                dv = dv_i;
        }
        t += dt;

        dt = dv_max/dv * dt;
        if (dt > dt_max)
            dt = dt_max;

        if (step % nStepsForReport == 0) {
            double px = 0, py = 0, pz = 0;
            double ek = 0;
            for (int i = 0 ; i < N; i++) {
                Particle *p = particles + i;
                px += p->vx;
                py += p->vy;
                pz += p->vz;
                ek += (p->vx*p->vx + p->vy*p->vy + p->vz*p->vz)/2;
            }


            for (unsigned i = 0; i < N; i++)
                u[i] = 0;
            cudaMemcpy(d_u, u, N*sizeof(double), cudaMemcpyHostToDevice);
            calcPotential<<<nBlocks, BLOCK_SIZE>>>(d_p, d_u, N);
            cudaMemcpy(u, d_u, N*sizeof(double), cudaMemcpyDeviceToHost);
            double ep = 0;
            for (unsigned i = 0; i < N; i++)
                ep += u[i];

            printf("%u %lf %lf %lf %lf %lf\n",
                step, t, dt, ep, ek, ek + ep);
        }

        if (step % nStepsForOutput == 0) {
            char fname[256];
            sprintf(fname, "out-%05u.txt", step);
            FILE *fout = fopen(fname, "w");
            for (unsigned i = 0; i < N; i++) {
                Particle *p = particles + i;
                fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
                        p->x, p->y, p->z, p->vx, p->vy, p->vz);
            }

            fclose(fout);
        }
    }

    printf("N=%d, Steps=%u Titer=%0.3lf s\n", N, nSteps, getTimer() / nSteps);
    free(particles);
    cudaFree(d_p);
    cudaFree(d_u);
}
