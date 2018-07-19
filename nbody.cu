#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

static clock_t t0;

void startTimer() {
    t0 = clock();
}

double getTimer() {
    return (clock() - t0) / (double)CLOCKS_PER_SEC;
}

#define BLOCK_SIZE 256

typedef struct {

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
} params_t;


static params_t params = {
    .n_steps = 1000*1000,
    .n_steps_for_report = 100,
    .n_steps_for_output = 1000,
    .dt_max = 0.01,
    .dv_max = 0.01,
    .dt_start = 0.001,
    .t_start = 0,
    .r2_eps = 0.01
};

typedef struct  {
    const char *name, *type;
    void *ptr;
} param_rec_t;

static const param_rec_t param_recs[] = {
    {   "n_steps",              "%u",   &params.n_steps             },
    {   "n_steps_for_report",   "%u",   &params.n_steps_for_report  },
    {   "n_steps_for_output",   "%u",   &params.n_steps_for_output  },
    {   "dt_max",               "%lf",  &params.dt_max              },
    {   "dv_max",               "%lf",  &params.dv_max              },
    {   "dt_start",             "%lf",  &params.dt_start            },
    {   "t_start",              "%lf",  &params.t_start             },
    {   "r2_eps",               "%lf",  &params.r2_eps              },
    {   NULL,   }
};

bool try_read_param(const char *line, const param_rec_t *rec) {
    const char
        *key = strstr(line, rec->name),
        *value = (key? strchr(key, '=') : NULL);
    return  value?  (sscanf(value + 1, rec->type, rec->ptr) == 1) : false;
}

void read_params() {

    FILE *fparam = fopen("params.txt", "r");
    if (!fparam)
        return;

    const unsigned max_line_len = 256;
    for(char buf[max_line_len] = {'\x0',}; fgets(buf, max_line_len, fparam);)
        for (const param_rec_t *rec = param_recs; rec->name; rec++)
            if (try_read_param(buf, rec))
                break;

    fclose(fparam);
}

void show_params(const bool to_stdout) {
    FILE *fo = to_stdout? stdout : stderr;
    for (const param_rec_t *rec = param_recs; rec->name; rec++) {
        (0 == strcmp("%lf", rec->type))?
            fprintf(fo, "%s = %lf\n", rec->name, *(double*)rec->ptr) :
        (0 == strcmp("%u", rec->type))?
            fprintf(fo, "%s = %u\n", rec->name, *(unsigned*)rec->ptr) :
            fprintf(fo, "%s = unknown\n", rec->name);
    }
}

typedef struct { double x, y, z, vx, vy, vz; } Particle;

__global__ void calcForces(Particle *p, double dt, unsigned N, double r2_eps) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double fx = 0, fy = 0, fz = 0;
        for (unsigned j = 0; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx*dx + dy*dy + dz*dz + r2_eps,
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

__global__
void calcPotential(Particle *p, double *u, unsigned N, double r2_eps) {
    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        double ui = 0;
        for (unsigned j = i + 1; j < N; j++) {
            const double
                dx = p[j].x - p[i].x,
                dy = p[j].y - p[i].y,
                dz = p[j].z - p[i].z,
                distSqr = dx*dx + dy*dy + dz*dz + r2_eps;
            ui += rsqrt(distSqr);
        }
        u[i] = ui;
    }
}


int main(const int argc, const char** argv) {

    if (argc >= 2 && strstr(argv[1], "-p"))
    {
        show_params(true);
        return 0;
    }

    read_params();
    show_params(false);

    Particle
        *particles = NULL,
        *old_particles = NULL;

    double *u = NULL;
    unsigned N = 0;

    {
        FILE *fin = fopen("input.txt", "r");
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

    double t = params.t_start, dt = params.dt_start;
    startTimer();
    for (unsigned step = 0; step < params.n_steps; step++) {

        cudaMemcpy(d_p, particles, N*sizeof(Particle), cudaMemcpyHostToDevice);
        calcForces<<<nBlocks, BLOCK_SIZE>>>(d_p, dt, N, params.r2_eps);
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

        dt = params.dv_max/dv * dt;
        if (dt > params.dt_max)
            dt = params.dt_max;

        if (step % params.n_steps_for_report == 0) {
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
            calcPotential<<<nBlocks, BLOCK_SIZE>>>(d_p, d_u, N, params.r2_eps);
            cudaMemcpy(u, d_u, N*sizeof(double), cudaMemcpyDeviceToHost);
            double ep = 0;
            for (unsigned i = 0; i < N; i++)
                ep += u[i];

            printf("%u %lf %lf %lf %lf %lf\n",
                step, t, dt, ep, ek, ek - ep);
        }

        if (step % params.n_steps_for_output == 0) {
            char fname[256];
            sprintf(fname, "out-%06u.txt", step);
            FILE *fout = fopen(fname, "w");
            for (unsigned i = 0; i < N; i++) {
                Particle *p = particles + i;
                fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
                        p->x, p->y, p->z, p->vx, p->vy, p->vz);
            }

            fclose(fout);
        }
    }

    printf("N=%d, Steps=%u Titer=%0.3lf s\n",
        N, params.n_steps, getTimer() / params.n_steps);
    free(particles);
    cudaFree(d_p);
    cudaFree(d_u);
}
