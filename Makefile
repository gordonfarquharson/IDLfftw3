CC = gcc
CINCLUDE = -I /usr/local/itt/idl/idl81/external/include
CFLAGS = -O2 -fPIC -Wall $(CINCLUDE)
LFLAGS = -lm -lfftw3f

CFILES = IDLfftw3.c
OFILES = IDLfftw3.o

TARGET = IDLfftw3.so

all:  	$(TARGET)

IDLfftw3.so: $(OFILES) 
	$(CC) -shared  $(CFLAGS) -o $(TARGET) $(OFILES) $(LFLAGS)

clean:
	rm -f $(OFILES) $(TARGET)

distclean:
	rm -f $(OFILES) $(TARGET) *~
