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


SRC = main.cpp readinputfiles.cpp processfile.cpp computerms.cpp ssa.cpp \
	  averagecycle.cpp smoothencoder.cpp fft.cpp scatterplot.cpp curvefit.cpp \
	  findroots.cpp printoutput.cpp spectrogram.cpp trackpeaks.cpp computemassflowrate.cpp \
	  computeballvalvearea.cpp bodeplot.cpp computeliquidfuelflowrate.cpp

OBJ = main.o readinputfiles.o processfile.o computerms.o ssa.o \
	  averagecycle.o smoothencoder.o fft.o scatterplot.o curvefit.o \
	  findroots.o printoutput.o spectrogram.o trackpeaks.o computemassflowrate.o \
	  computeballvalvearea.o bodeplot.o computeliquidfuelflowrate.o

processPC: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o processPC $(LDLIBS)
main.o: main.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c main.cpp -o $@ $(LDLIBS)
readinputfiles.o: readinputfiles.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c readinputfiles.cpp -o $@ $(LDLIBS)
processfile.o: processfile.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c processfile.cpp -o $@ $(LDLIBS)
ssa.o: ssa.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c ssa.cpp -o $@ $(LDLIBS)
computerms.o: computerms.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c computerms.cpp -o $@ $(LDLIBS)
averagecycle.o: averagecycle.cpp findroots.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c averagecycle.cpp -o $@ $(LDLIBS)
fft.o: fft.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c fft.cpp -o $@ $(LDLIBS)
scatterplot.o: scatterplot.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c scatterplot.cpp -o $@ $(LDLIBS)
smoothencoder.o: smoothencoder.cpp ssa.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c smoothencoder.cpp -o $@ $(LDLIBS)
curvefit.o: curvefit.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c curvefit.cpp -o $@ $(LDLIBS)
findroots.o: findroots.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c findroots.cpp -o $@ $(LDLIBS)
printoutput.o: printoutput.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c printoutput.cpp -o $@ $(LDLIBS)
spectrogram.o: spectrogram.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c spectrogram.cpp -o $@ $(LDLIBS)
trackpeaks.o: trackpeaks.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c trackpeaks.cpp -o $@ $(LDLIBS)
computemassflowrate.o: computemassflowrate.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c computemassflowrate.cpp -o $@ $(LDLIBS)
computeballvalvearea.o: computeballvalvearea.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c computeballvalvearea.cpp -o $@ $(LDLIBS)
bodeplot.o: bodeplot.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c bodeplot.cpp -o $@ $(LDLIBS)
computeliquidfuelflowrate.o: computeliquidfuelflowrate.cpp 
	$(CC) $(CFLAGS) $(LDFLAGS) -c computeliquidfuelflowrate.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm $(OBJ)
	rm processPC
