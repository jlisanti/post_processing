# Compiler options

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CC := g++
endif

ifeq ($(UNAME), Linux)
CC := g++-4.8
endif

LDFLAGS := -L/usr/local/lib \
           -I/usr/local/include \
           -I/Users/jlisanti/Documents/program_sources/boost_1_67_0 \
		   -I/Users/mohammadqatomah/boost_1_67_0 \
		   -L/opt/local/lib \
           -I/opt/local/include \


LDLIBS := -ldatafileinput -lgsl -lgslcblas -lm -lmenu
CFLAGS := -Wall -DHAVE_INLINE -v -std=c++11

SRC = main.cpp

OBJ = main.o

valveFunction: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o valveFunction $(LDLIBS)
main.o: main.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c main.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm $(OBJ)
	rm valveFunction 
