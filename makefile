CC=gcc
CFLAGS=-I. -Wall
LIBS=-lm

DEPS = libjonswap.h
OBJ = jonswap.o libjonswap.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

jonswap: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

