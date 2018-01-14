# Makefile for Selector

# Compiler
CC = g++

# Compiler options
CFLAGS = -g -Wall -pedantic -fopenmp -std=c++0x

# all object files
SEL_OBJECTS = selector_user.o selector.o selector_internal.o

selector : $(SEL_OBJECTS)
	$(CC) $(CFLAGS) $(SEL_OBJECTS) -o lex -lm

selector_internal.o : selector_internal.c selector_internal.h selector.h selector_user.h
	$(CC) $(CFLAGS) -c selector_internal.c 

selector_user.o : selector_user.c selector_user.h selector.h
	$(CC) $(CFLAGS) -c selector_user.c

selector.o : selector.c selector.h selector_user.h selector_internal.h
	$(CC) $(CFLAGS) -c selector.c

clean:
	rm -f *~ *.o
