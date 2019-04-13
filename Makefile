CFLAGS = -g -Wall -Wno-deprecated

all: clean main

main: image.h rgb.h rgb.o image.o ReadImage.o ReadImageHeader.o WriteImage.o \
	main.cpp
	g++ -o main $(CFLAGS) image.o rgb.o ReadImage.o ReadImageHeader.o \
					WriteImage.o main.cpp

ReadImage.o:	image.h rgb.h ReadImage.cpp
	g++ -c $(CFLAGS) ReadImage.cpp

ReadImageHeader.o:	image.h rgb.h ReadImageHeader.cpp
	g++ -c $(CFLAGS) ReadImageHeader.cpp

WriteImage.o:	image.h rgb.h WriteImage.cpp
	g++ -c $(CFLAGS) WriteImage.cpp

image.o:	image.h rgb.h image.cpp
	g++ -c $(CFLAGS) image.cpp

rgb.o:	rgb.h rgb.cpp
	g++ -c $(CFLAGS) rgb.cpp

clean:
	rm -f main *.o