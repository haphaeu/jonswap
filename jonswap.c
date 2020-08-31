#include "libjonswap.h"

int main(int argc, char *argv[]) {

    double hs, tp;
    double w1, w2;
    double to, tf, ts;
    double *t, *w, *pm, *js, *amp, *phi, *tt, *eta;
    double (*fgamma)(double, double);
    int nharms, tt_size, seed;
    short show_spectrum, show_timetrace;
    
    // hard coded inits
    w1 = 0.2;
    w2 = 6.28;
    to = 0.0;

    if (--argc < 2) {
        printf("Error: invalid arguments. Use:\n\n");
        printf("  jonswap hs tp [-n num_harmonics] [-g gamma] [-s seed] [-d duration] [-t timestep] [-v] [-h]\n");
        printf("\n  hs: significant wave height in meters.");
	printf("\n  tp: spectral peak period in seconds");
        printf("\n  -n: number of harmonics used to discretise the spectrum.");
        printf("\n  -g: value for gamma. If ommited, DNV is used.");
        printf("\n  -s: seed number for phase randomisation.");
	printf("\n  -d: duration - timetrace will be shown.");
	printf("\n  -t: time step. default is 0.1 seconds.");	
        printf("\n  -v: verbose output");
        printf("\n  -h: hide spectrum");
        printf("\n");
        exit(1);
    }
    hs = atof(*++argv); argc--;
    tp = atof(*++argv); argc--;
    
    nharms = 200;
    fgamma = jonswap_gamma;
    seed = SEED;
    tf = 0.0;
    ts = 0.1;
    show_timetrace = 0;
    show_spectrum = 1;
    char c;
    while (argc-- > 0 && (*++argv)[0] == '-') {
        c = *++argv[0];
        switch (c) {
        case 'n':
            if (!argc-- || (nharms = atoi(*++argv)) == 0) {
            printf("Error: invalid nharms argument.\n");
            exit(2);
            }
            break;
        case 'g':
            if (!argc-- || (custom_gamma_val = atof(*++argv)) == 0.0) {
            printf("Error: invalid gamma argument.\n");
            exit(3);
            }
            fgamma = custom_gamma_func;     
            break;
        case 's':
            if (!argc-- || (seed = atoi(*++argv)) == 0) {
            printf("Error: invalid seed argument.\n");
            exit(4);
            }
            break;
	case 'd':
	    if (!argc-- || (tf = atof(*++argv)) == 0.0) {
            printf("Error: invalid duration argument.\n");
            exit(5);
            }
	    show_timetrace = 1;
            break;
	case 't':
	    if (!argc-- || (ts = atof(*++argv)) == 0.0) {
            printf("Error: invalid time step argument.\n");
            exit(6);
            }
	    show_timetrace = 1;
            break;
        case 'v':
            verbose = 1;
            break;
        case 'h':
            show_spectrum = 0;
            break;
        default:
            printf("Warning: illegal input option -%c.\n", c);
            break;
        }
    }
    if (++argc > 0)
        printf("Warning: %d unprocessed input argument%s.\n", argc, argc == 1 ? "" : "s");

    if (verbose) printf("Hs %.2f m\nTp %.2f s \nDuration %.2f s \nStep %.2f s \nnharms %d\nseed %d\n", hs, tp, tf, ts, nharms, seed);


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
        printf("\nSpectrum\n\n");
        printf("%10s %10s %12s %12s %12s %10s \n", "T", "w", "PM", "JS", "amp", "phi");
        printf("%10s %10s %12s %12s %12s %10s \n", "[s]", "[rd/s]", "[m2s/rd]", "[m2s/rd]", "[m]", "[rd]");
        for (int i = 0; i < nharms; ++i)
            printf("%10.3f %10.3f %12.5f %12.5f %12.5f %10.3f\n", t[i], w[i], pm[i], js[i], amp[i], phi[i]);
    }
    if (show_timetrace) {
        printf("\nTime History\n\n");
        printf("%10s %10s \n", "Time", "Elevation");
        for (int j = 0; j < tt_size; ++j)
            printf("%10.2f %10.3f\n", tt[j], eta[j]);
    }
    
    return 0;
}

