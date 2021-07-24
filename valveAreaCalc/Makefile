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
           -I/Users/Joel/Documents/programSourceFiles/boost_1_72_0 \
           -I/Users/joellisanti/Documents/programSourceFiles/boost_1_72_0 \
		   -I/Users/mohammadqatomah/boost_1_72_0 \
		   -L/opt/local/lib \
           -I/opt/local/include \


LDLIBS := -ldatafileinput -lgsl -lgslcblas -lm -lmenu
CFLAGS := -Wall -DHAVE_INLINE -v -std=c++11


SRC = main.cpp computeballvalvearea.cpp

OBJ = main.o computeballvalvearea.o

valveCalc: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o valveCalc $(LDLIBS)
main.o: main.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c main.cpp -o $@ $(LDLIBS)
computeballvalvearea.o: computeballvalvearea.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c computeballvalvearea.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm $(OBJ)
	rm valveCalc
