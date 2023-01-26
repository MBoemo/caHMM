CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS =
CXXFLAGS = -fopenmp -std=c++11 $(DEBUG)
CFLAGS = -std=c99 -O2 $(DEBUG)

MAIN_LIBRARY = libcaHMM.a

all: depend $(MAIN_LIBRARY)

SUBDIRS = src
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))

#generate object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

depend: .depend

.depend: $(CPP_SRC) $(C_SRC)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

#compile each object
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(H5_INCLUDE) -fPIC $<

#compile the main executable
$(MAIN_LIBRARY): $(CPP_OBJ) $(C_OBJ)
	ar rcs -o $@ $(CPP_OBJ) $(C_OBJ) $(LIBFLAGS)

clean:
	rm -f $(MAIN_LIBRARY) $(CPP_OBJ) $(C_OBJ)
