compiler = gcc
cflags = -std=c99 -lm

headers = $(wildcard *.h)
sources = $(wildcard *.c)
objects = $(sources:.c=.o)

executables: arnoldi

%.o: %.c $(headers)
	$(compiler) -c -o $@ $< $(cflags)

arnoldi: $(objects)
	$(compiler) -o $@ $^ $(cflags)

clean:
	rm -f arnoldi *.o
