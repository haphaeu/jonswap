/*

* TODO
 - [X] Re-think PERIOD_START and PERIOD_END constants
 - [X] Re-define spectrum based on frequency and sort out a-/de-scending order of omega array.
 - [X] Implement mandatory input arguments for
   - [X] simulation duration
   - [X] time step 
 - [X] Implement a function to parse input arguments to handle optional input arguments
   - [X] number of harmonics
   - [X] gamma parameter
   - [X] random seed

 - Think batch mode - maybe convert data into structs

c:/Users/rarossi/portable/Cpp/App/Dev-Cpp/MinGW64/bin/gcc -std=c11 jonswap.c -o jonswap

gcc -lm jonswap.c -o jonswap.out

* JONSWAP spectrum calculations

*/

#include "libjonswap.h"

short verbose = 0;

double jonswap_gamma(double hs, double tp) {
    /* 
    Return gamma factor for a JONSWAP spectrum according to DNV
    */
    double chk;

    chk = tp / sqrt(hs);

    if (chk < 3.6)
    return 5.0;
    else if (chk > 5.0)
    return 1.0;
    else
    return exp(5.75 - 1.15*chk);
}

double custom_gamma_val = 0.0;
double custom_gamma_func(double hs, double tp) {
    return custom_gamma_val;
}

double* freq_domain(double w1, double w2, int length) {
    /*
    Returns an array of angular frequencies [rd/s] from w1 to w2, of size length.
     */
    double *f = (double *) malloc(length * sizeof(double));

#if FREQUENCY_DISCRETISATION == ARITHMETIC 
    // Arithmetic progression on the frequency f[i+i] = f[i] + step
    // Leads to poor discretisation
    double step = (w2 - w1) / (length - 1);
    for (int i = 0; i < length; ++i)
        *(f+i) = w1 + i*step;
#elif FREQUENCY_DISCRETISATION == GEOMETRIC
    // Geometric progression f[i+1] = r * f[i]
    // Better discretisation in the region of interest.
    double r = pow(w2/w1, 1.0/(length-1));
    f[0] = w1;
    for (int i = 1; i < length; ++i)
        *(f+i) = *(f+i-1) * r;
#else // error
    free(f);
    return NULL;
#endif
    return f;
}

double* period_domain(double *w, int length) {
    /*
    Returns an array of periods [s] corresponding to the angular frequencies from the array w.
     */
    double *t = (double*) malloc(length * sizeof(double));
    for (int i = 0; i < length; ++i)
        *(t+i) = TWOPI / *(w+i);
    return t;
}

double* pierson_moskowitz_spectrum(double hs, double tp, int length, double *w) {
    /*
    Return Pierson Moskowitz spectrum.

    hs     : significant wave height
    tp     : peak period
    length : length of the array
    w      : angular frequencies, or NULL to calculate it.
     */
    double wp, pm_cte, *s;

    wp = TWOPI / tp;
    pm_cte = 0.3125 * pow(hs, 2) * pow(wp, 4);

    if (verbose) printf("Pierson Moskowitz\n  wp %.2f\n  cte %.2f\n", wp, pm_cte);

    s = (double*) malloc(length * sizeof(double));

    for (int i = 0; i < length; ++i)
        *(s+i) = pm_cte * pow(*(w+i), -5) * exp(-1.25 * pow(*(w+i)/wp, -4));

    return s;
}

double* jonswap_spectrum(double hs, double tp, double (*fgamma)(double, double),
             int length, double *w, double *pm) {
    /*
    Return JONSWAP spectrum.

    hs     : significant wave height
    tp     : peak period
    fgamma : function(hs, tp) retuning gamma parameter
    length : length of the array
    w      : angular frequencies
    pm     : Pierson Moskowitz spectrum

    */
    double wp, gamma, norm, *s;
    wp = TWOPI / tp;
    gamma = fgamma(hs, tp);
    norm = 1.0 - 0.287 * log(gamma);

    if (verbose) printf("JONSWAP\n  wp %.2f\n  gamma %.2f\n  norm %.2f\n", wp, gamma, norm);

    s = (double*) malloc(length * sizeof(double));

    for (int i = 0; i < length; ++i)  // I might regret this later, but I'm using
        *(s+i) = norm * *(pm+i) *     // pointer shift for arrays instead of indexes, i.e:
            pow(gamma,                // *(pm+i) == pm[i]
            exp(-0.5 * pow((*(w+i) - wp)/wp/ (*(w+i) < wp ? 0.07 : 0.09), 2)));
    return s;
}

double spectral_moment(int n, double *s, double *w, int length) {
    /*
    Return the n-th spectral moment.

    n      : order of the spectral moment to be calculated
    s      : spectrum array
    w      : angular frequency array
    length : length of the s and w arrays
     */
    double m;
    for (int i = 0; i < length - 1; ++i)
        m += pow((*(w+i+1) + *(w+i))/2.0, n) * (*(w+i+1) - *(w+i)) * (*(s+i) + *(s+i+1));/* / 2.0; */
    return m/2.0;
}

double* jonswap_component_amplitude(double *s, double *w, int length) {
    /*
    Return the amplitude of a sinusoidal signal for each frequency component
    of a discretizes spectrum.
    
    s      : spectrum array
    w      : angular frequency array
    length : length of the arrays

    Ref: https://ceprofs.civil.tamu.edu/jzhang/ocen300/statistics-spectrum.pdf

    */
    double *amp = (double*) malloc(length * sizeof(double));
    for (int i = 0; i < length - 1; ++i)
        *(amp+i) = sqrt(2 * *(s+i) * (*(w+i+1) - *(w+i)));
    *(amp+length-1) = 0.0;
    return amp;
}

double* phases(int length, int seed) {
    double *phi = (double*) malloc(length  * sizeof(double));
    srand(seed);
    for (int i = 0; i < length; ++i)
        // random double between -pi and pi
        *(phi+i) = ((double)rand()/RAND_MAX*2.0 - 1.0) * PI;
    return phi;
}
double* time_domain(double to, double tf, double ts, int *length) {
    /*
    Return an array of times domain [s] to be used to calculate the wave elevation.

    Note that since input is time-step, the array length must be passed
    as a pointer and it is where the length will be stored. 
    */
    double *td;
    *length = (int) ((tf - to) / ts);
    td = (double*) malloc(*length * sizeof(double));
    for (int i = 0; i < *length; ++i)
        *(td+i) = to + i*ts;
    return td;    
}
double *wave_elevation(double *amp, double *w, double *phi, double *t, int len_spectrum, int len_tt) {
    double *eta;
    eta = (double*) malloc(len_tt * sizeof(double));
    for (int i = 0; i < len_tt; ++i) {
    eta[i] = 0.0;
    for (int j = 0; j < len_spectrum; ++j)
        // ok, here I got tired of pointer notation =PPP
        //eta[i] += amp[j] * cos(w[j]*t[i] - phi[j]);
        *(eta+i) += *(amp+j) * cos(*(w+j) * *(t+i) - *(phi+j));
    }
    return eta;
}
