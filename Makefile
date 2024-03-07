#
# Makefile to compile LFmap
#
#
#
CC        = gcc

CFLAGS = -I/usr/local/include -I/home/emilp/Desktop/work/cfitsio
LDFLAGS  = -L/usr/local/lib -L/home/emilp/Desktop/work/cfitsio
LDLIBS   =  -lcfitsio -lX11 -lm

LFmap: LFmap.o UtilLFmap.o
	cc $(LDFLAGS) LFmap.o UtilLFmap.o -o LFmap $(LDLIBS)

LFmap.o: UtilLFmap.h

UtilLFmap.o: UtilLFmap.h


clean:
	rm -f *.o core* LFmap
