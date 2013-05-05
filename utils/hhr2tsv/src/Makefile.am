bin_PROGRAMS=hhr2tsv
hhr2tsv_SOURCES=convert.h convert.c main.c 

convert.o:
	$(CC) -DEXIT_CODE=0 -c convert.c

main.o: convert.o
	$(CC) -DEXIT_CODE=0 -c main.c
