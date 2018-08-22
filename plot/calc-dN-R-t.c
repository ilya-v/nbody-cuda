#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

    if (argc < 2) return 1;
    FILE *fin = fopen(argv[1], "r");
    if (!fin) return 2;

    double xc = 0.0, yc = 0.0, zc = 0.0;
    unsigned N = 0;

    for(char line[2048]; fgets(line, sizeof(line), fin);) {
        double x, y, z;
        if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
            return 3;
        xc += x; yc += y; zc +=z; N++;
    }
    xc /= N; yc /= N; zc /= N;

    const double r_max = 1000.0, dr = 1;
    const unsigned nr = (int)floor(r_max/dr) + 2;
    unsigned *dN = calloc(nr, sizeof(unsigned));

    rewind(fin);

    for(char line[2048]; fgets(line, sizeof(line), fin);) {
        double x,y,z;
        if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
            return 3;
        x -= xc; y -= yc; z -= zc;
        double r = sqrt(x*x+y*y+x*x);
        unsigned int i = (int)floor(r/dr);
        i = ( (i < nr - 1)? i : nr -1);
        dN[i] += 1;
    }

    for (unsigned i = 0; i < nr; i++)
        printf("%u\n", dN[i]);
    return 0;
}
