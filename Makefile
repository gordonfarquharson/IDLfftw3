CC = gcc
CINCLUDE = -I /usr/local/itt/idl/idl81/external/include
CFLAGS = -O2 -fPIC -Wall $(CINCLUDE)
LFLAGS = -lm -lfftw3f

CFILES = IDLfftw3.c
OFILES = IDLfftw3.o

all:  	libIDLfftw.so

libIDLfftw.so: $(OFILES) 
	$(CC) -shared  $(CFLAGS) -o libIDLfftw.so $(OFILES) $(LFLAGS)

clean:
	rm -f *.o libIDLfftw.so

distclean:
	rm -f *.o *~ libIDLfftw.so