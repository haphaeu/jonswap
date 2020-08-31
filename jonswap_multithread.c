/*
Mock-up implementation using multithreading to benchmark performance.

This does not do anything specific rather than calculating various
spectra and timetraces with very refined steps.

Compile command:
gcc jonswap_multithread.c -o jonswap_multithread -lm -lpthread -Wno-implicit-function-declaration

#+AUTHOR Rafael Rossi
#+DATE 25-Oct-2019

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <pthread.h>
#include <time.h>

#include "libjonswap.h"

#define MAX_NUM_THREADS 99

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
	if (verbose) {
	    printf("\nSpectrum\n\n");
	    printf("%10s %10s %12s %12s %12s %10s \n", "T", "w", "PM", "JS", "amp", "phi");
	    printf("%10s %10s %12s %12s %12s %10s \n", "[s]", "[rd/s]", "[m2s/rd]", "[m2s/rd]", "[m]", "[rd]");
	    for (int i = 0; i < nharms; ++i)
		printf("%10.3f %10.3f %12.5f %12.5f %12.5f %10.3f\n", t[i], w[i], pm[i], js[i], amp[i], phi[i]);
	}
    }
    if (show_timetrace) {
	if (verbose) {
	    printf("\nTime History\n\n");
	    printf("%10s %10s \n", "Time", "Elevation");
	    for (int j = 0; j < tt_size; ++j)
		printf("%10.2f %10.3f\n", tt[j], eta[j]);
	}
    }
    
    printf("Done with Hs %.1f Tp %4.1f dt %.2f seed %ld\n", hs, tp, ts, seed);
    
    return NULL;
}

