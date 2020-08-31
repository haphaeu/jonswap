/*
Calculates a discretised JONSWAP spectrum based on input Hs and Tp.

Also calculates a time history of wave elevation according to that spectrum.

#+AUTHOR Rafael Rossi
#+DATE 06-Jun-2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI            3.14159265358979323846
#define TWOPI         6.28318530717958647692
#define SEED          12547

#define ARITHMETIC                  01 
#define GEOMETRIC                   02
#define FREQUENCY_DISCRETISATION    GEOMETRIC

short verbose;

double jonswap_gamma(double hs, double tp);
double custom_gamma_val;            // dummy variable and function to handle user custom...
double custom_gamma_func(double, double); // ...value for the gamma parameter
double* freq_domain(double w1, double w2, int length);
double* period_domain(double *w, int length);
double* pierson_moskowitz_spectrum(double hs, double tp, int length, double *w);
double* jonswap_spectrum(double hs, double tp, double (*fgamma)(double, double),
             int length, double *w, double *pm);
double spectral_moment(int n, double *s, double *w, int length);
double* jonswap_component_amplitude(double *s, double *w, int length);
double* phases(int length, int seed);
double* time_domain(double to, double tf, double ts, int *length);
double* wave_elevation(double *amp, double *w, double *phi, double *t, int len_spectrum, int len_tt);
