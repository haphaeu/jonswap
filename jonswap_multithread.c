/*
Mock-up implementation using multithreading to benchmark performance.

This does do anything specific rather than calculating various
spectra and timetraces with very refined steps.


gcc jonswap_multithread.c gnuplot_i.c -o jonswap_multithread -lm -lpthread -Wno-implicit-function-declaration

#+AUTHOR Rafael Rossi
#+DATE 25-Oct-2019

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "gnuplot_i.h"

#include <pthread.h>
#include <time.h>

#define PI            3.14159265358979323846
#define TWOPI         6.28318530717958647692
#define SEED          12547

#define ARITHMETIC                  01 
#define GEOMETRIC                   02
#define FREQUENCY_DISCRETISATION    GEOMETRIC

#define MAX_NUM_THREADS 99

short verbose = 0;

double jonswap_gamma(double hs, double tp);
double custom_gamma_val = 0.0;            // dummy variable and function to handle user custom...
double custom_gamma_func(double, double); // ...value for the gamma parameter
double* freq_domain(double w1, double w2, int length);
double* period_domain(double *w, int length);
double* pierson_moskowitz_spectrum(double hs, double tp, int length, double *w);
double* jonswap_spectrum(double hs, double tp, double (*fgamma)(double, double),
             int length, double *w, double *pm);
double spectral_moment(int n, double *s, double *w, int length);
double* jonswap_component_amplitude(double *s, double *w, int length);
double* phases(int length, long int seed);
double* time_domain(double to, double tf, double ts, int *length);
double* wave_elevation(double *amp, double *w, double *phi, double *t, int len_spectrum, int len_tt);

void plot(double *x, double *y, int npoints);

void* worker_thread(void *vargp);
double randbetween(double a, double b);


struct WorkerThreadArgs {
    double hs;
    double tp;
    double tf;
    double ts;
    int seed;
} thread_args[MAX_NUM_THREADS]; 


int main(int argc, char *argv[]) {
    
    pthread_t thread_id[MAX_NUM_THREADS];
    int num_threads;

    struct timespec start, finish;
    double elapsed;

    srandom(time(NULL));

    if(argc == 1) {
	// default option
	num_threads = 3;
    } else {
	int n = atoi(argv[1]);
	num_threads = n < MAX_NUM_THREADS ? n : MAX_NUM_THREADS;
    }									    

    // define a queue of random args for the threads
    for(int i=0; i<num_threads; ++i) {
	thread_args[i] = (struct WorkerThreadArgs) {
	    .hs=randbetween(1.5, 3.5),
	    .tp=randbetween(6, 14),
	    .tf=10000,
	    .ts=randbetween(0.01, 0.25),
	    .seed=random()
	};
    }

    ////////////////////////////////////////////////////////////////////////
    /////////// Running multi threaded /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    printf("Running %d thread%c.\n", num_threads, num_threads > 1 ? 's':' ');
    clock_gettime(CLOCK_MONOTONIC, &start);

    // create and launch all threads
    for(int i=0; i<num_threads; ++i)
	pthread_create(&thread_id[i], NULL, worker_thread, thread_args+i);

    // wait until all threads are done
    for(int i=0; i<num_threads; ++i)    
        pthread_join(thread_id[i], NULL);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Done with all threads - %.3fs\n", elapsed);

    ////////////////////////////////////////////////////////////////////////
    /////////// Running single threaded/////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    printf("Running sequentially.\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    // create and launch all threads
    for(int i=0; i<num_threads; ++i) {
	pthread_create(&thread_id[i], NULL, worker_thread, thread_args+i);
        pthread_join(thread_id[i], NULL); // join after create will wait until it is done
    }
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Done sequentially - %.3fs\n", elapsed);

    //    printf("Press Enter to exit.\n");
    //    char input[5];
    //    fgets(input, sizeof input, stdin);
    exit(0); 

}


double randbetween(double a, double b) {
    return a + (double) rand() /(double) (RAND_MAX/(b - a));
}


void* worker_thread(void *vargp) {

    struct WorkerThreadArgs* args = (struct WorkerThreadArgs *) vargp; 

    double hs = args->hs;
    double tp = args->tp;
    double tf = args->tf;
    double ts = args->ts;
    long int seed = args->seed;
    double to;
    double w1, w2;
    double *t, *w, *pm, *js, *amp, *phi, *tt, *eta;
    double (*fgamma)(double, double);
    int nharms, tt_size;
    short show_spectrum, show_timetrace;
    
    // handler for the gnuplot interface
    gnuplot_ctrl *h1, *h2;
    h1 = gnuplot_init(); h2 = gnuplot_init();
    gnuplot_setstyle(h1, "lines"); gnuplot_setstyle(h2, "lines");

    // hard coded inits
    w1 = 0.2;
    w2 = 6.28;
    to = 0.0;

    nharms = 200;
    fgamma = jonswap_gamma;
    show_timetrace = 0;
    show_spectrum = 0;
  
    if (verbose) printf("Hs %.2f m\nTp %.2f s \nDuration %.2f s \nStep %.2f s \nnharms %d\nseed %ld\n",
			hs, tp, tf, ts, nharms, seed);

    w = freq_domain(w1, w2, nharms);
    t = period_domain(w, nharms);
    pm = pierson_moskowitz_spectrum(hs, tp, nharms, w);
    js = jonswap_spectrum(hs, tp, fgamma, nharms, w, pm);
    amp = jonswap_component_amplitude(js, w, nharms);
    phi = phases(nharms, seed);
    tt = time_domain(to, tf, ts, &tt_size);
    eta = wave_elevation(amp, w, phi, tt, nharms, tt_size);

    if (verbose)
        // checking using the relationship zero-th moment is Hs squared divided by 16
        // m0 = Hs**2 / 16 ==> Hs = 4*sqrt(m0)
        printf("Check Hs %.2f m\n", 4*sqrt(spectral_moment(0, js, w, nharms)));
    if (show_spectrum){
        gnuplot_plot_xy(h1, t, js, nharms, "Spectrum");
	if (verbose) {
	    printf("\nSpectrum\n\n");
	    printf("%10s %10s %12s %12s %12s %10s \n", "T", "w", "PM", "JS", "amp", "phi");
	    printf("%10s %10s %12s %12s %12s %10s \n", "[s]", "[rd/s]", "[m2s/rd]", "[m2s/rd]", "[m]", "[rd]");
	    for (int i = 0; i < nharms; ++i)
		printf("%10.3f %10.3f %12.5f %12.5f %12.5f %10.3f\n", t[i], w[i], pm[i], js[i], amp[i], phi[i]);
	}
    }
    if (show_timetrace) {
        gnuplot_plot_xy(h2, tt, eta, tt_size, "Time History");
	if (verbose) {
	    printf("\nTime History\n\n");
	    printf("%10s %10s \n", "Time", "Elevation");
	    for (int j = 0; j < tt_size; ++j)
		printf("%10.2f %10.3f\n", tt[j], eta[j]);
	}
    }
    
    //    printf("Press Enter to exit.\n");
    //    char input[5];
    //    fgets(input, sizeof input, stdin);
    //    gnuplot_close(h1); gnuplot_close(h2);

    printf("Done with Hs %.1f Tp %4.1f dt %.2f seed %ld\n", hs, tp, ts, seed);
    
    return NULL;
}

double jonswap_gamma(double hs, double tp) {
    /* 
    Return gamma factor for a JONSWAP spectrum according to DNV
    */
    double chk, gamma;

    chk = tp / sqrt(hs);

    if (chk < 3.6)
    return 5.0;
    else if (chk > 5.0)
    return 1.0;
    else
    return exp(5.75 - 1.15*chk);
}

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
        m += pow((*(w+i+1) - *(w+i))/2.0, n) * (*(w+i+1) - *(w+i)) * (*(s+i) + *(s+i+1));/* / 2.0; */
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

double* phases(int length, long int seed) {
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

