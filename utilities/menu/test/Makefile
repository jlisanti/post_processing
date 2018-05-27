# Compiler options

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CC := g++
endif

ifeq ($(UNAME), Linux)
CC := g++-4.8
endif

LDFLAGS := -L/usr/local/lib -I/usr/local/include
LDLIBS := -lmenu
CFLAGS := -Wall -DHAVE_INLINE -v -std=c++11


SRC = main.cpp

OBJ = main.o

test-menu: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o test-menu $(LDLIBS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c main.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm $(OBJ)
	rm test-menu
