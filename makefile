CC=gcc
CFLAGS=-I. -Wall
LIBS=-lm

DEPS = libjonswap.h
OBJ = jonswap.o libjonswap.o
OBJ2 = plot.o libjonswap.o gnuplot_i.o
OBJ3 = jonswap_multithread.o libjonswap.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

jonswap: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

plot: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) -Wno-implicit-function-declaration

jonswap_multithread: $(OBJ3)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) -lpthread
