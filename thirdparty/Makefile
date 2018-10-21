FC	=	gfortran
CC	=	g++
INC	=	/usr/include
LIB	=	/usr/lib64
RFM	= rfm-4.33
DISORT	= cdisort-2.1.3a

all: hitbin kcoeff mppnccombine
	make -C rfm-4.33 
	make -C cdisort-2.1.3a -f Makefile_cdisort
hitbin: hitbin.f
	$(FC) hitbin.f -o hitbin
kcoeff: kcoeff.cpp
	$(CC) -O3 kcoeff.cpp -I $(INC) -L$(LIB) -lnetcdf -o kcoeff
mppnccombine: mppnccombine.c
	gcc -O3 mppnccombine.c -I $(INC) -L$(LIB) -lnetcdf -o mppnccombine
clean:
	make -C rfm-4.33 clean
	make -C cdisort-2.3.1a -f Makefile_cdisort clean
	rm -f hitbin kcoeff mppnccombine
.PHONY:
	all clean
