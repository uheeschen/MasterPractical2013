hhr2tsvdir=../
hhr2tsv_PROGRAMS=hhr2tsv
hhr2tsv_SOURCES=convert.c main.c
hhr2tsv_CFLAGS=-I../include

convert.o:
	$(COMPILE) -DEXIT_CODE=0 -c convert.c

main.o: convert.o
	$(COMPILE) -DEXIT_CODE=0 -c main.c
