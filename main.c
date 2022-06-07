#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define g 9.81f
#define m 1.0f
#define alpha 0.5f

void pochodne(double t, double *s, double *k) {
    k[0]=s[2];
    k[1]=s[3];
    k[2]=-g*pow(cos(alpha), 2)*sin(s[0]) / (sin(alpha)*s[1]) - 2*s[2]*s[3]/s[1];
    k[3]=pow(sin(alpha), 2)*s[1]*pow(s[2], 2) - g*sin(alpha)*pow(cos(alpha), 2)*(1-cos(s[0]));
}

void rk4_vec(double t, double dt, int n, double *s, void (*f)(double , double *, double  *)) {
    #define M 1000
    static double k1[M], k2[M], k3[M], k4[M], w[M];
    int i;
    // Kopia tablicy
    for(i = 0; i < n; i++) {
        w[i] = s[i];
    }
    f(t, w, k1);
    // get rk1
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k1[i];
    }
    f(t+dt/2, w, k2);
    // get rk2
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k2[i];
    }
    f(t+dt/2, w, k3);
    // get rk3
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt*k3[i];
    }
    f(t+dt, w, k4);
    // get rk4, save results
    for(i = 0; i < n; i++) {
        s[i]= s[i] + dt/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

void changeSystem(double* s, double* s2) {
    double rho = s[1]*tan(alpha);
    double x = rho*cos(s[0]);
    double y = rho*sin(s[0]);
    double z = s[1];

    double theta = 3.14/2 - alpha;
    s2[0] = cos(theta)*x + sin(theta)*z;
    s2[1] = y;
    s2[2] = -sin(theta)*x + cos(theta)*z;
}

int main() {
    // inicjalizacja  parametrów:
    int n = 4;              // ilość zmiennych w układzie RRZ1
    double dt = 0.1;
    double tmax = 50;
    int N = (int)tmax/dt;    // ilość kroków  czasowych
    double t = 0;

    void (*f)(double, double *, double *);
    f = pochodne;

    double *s;
    s = (double *)malloc(n*sizeof(double));    // tablica  rozwiązań

    double *s2;
    s2 = (double *)malloc(3*sizeof(double));

    // warunki początkowe
    s[0] = 1.1; // radians
    s[1] = 1.0;
    s[2] = g;
    s[3] = g;

    // plik
    FILE *fp;
    fp = fopen("wyniki.txt", "w");
    if(fp == NULL) {
        printf("Error: No file found\n");
        return EXIT_FAILURE;
    }

    // transformacja do układu laboratoryjnego
    FILE *fp2;
    fp2 = fopen("laboratoryjny.dat", "w");
    if(fp == NULL) {
        printf("Error: No file found\n");
        return EXIT_FAILURE;
    }

    // pętla symulacji
    for(size_t i = 1; i <= N; i++) {
        rk4_vec(t, dt, n, s, f);

        double fi = s[0];
        double z = s[1];
        double dfi = s[2];
        double dz = s[3];
        double e = (1/2.0f)*(pow(tan(alpha), 2) * pow(z, 2) * pow(dfi, 2) + pow(dz, 2)/pow(cos(alpha), 2)) + g*z*sin(alpha)*(1-cos(fi));
        t=t+dt;
        fprintf(fp, "%f %f %f %f %f %f\n", t, s[0], s[1], s[2], s[3], e);
        // trajektoria w układzie obserwatora
        changeSystem(s, s2);
        fprintf(fp2, "%f %f %f\n", s2[0], s2[1], s2[2]);
    }
    fclose(fp);
    fclose(fp2);
    return EXIT_SUCCESS;
}
