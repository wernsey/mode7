# Makefile for Linux and Windows (MinGW)
CC=gcc
CFLAGS=-c -Wall -Isdl -I.
LDFLAGS=-lm

SOURCES=sdl/pocadv.c bmp.c render.c mode7.c obj.c
OBJECTS=$(SOURCES:.c=.o)

# Detect operating system:
# More info: http://stackoverflow.com/q/714100
ifeq ($(OS),Windows_NT)
  EXECUTABLE=mode7.exe
  LDFLAGS+=-mwindows
else
  EXECUTABLE=mode7
endif

ifeq ($(BUILD),debug)
  # Debug
  CFLAGS += -O0 -g -I/local/include
  LDFLAGS +=
else
  # Release mode
  CFLAGS += -O2 -DNDEBUG -I/local/include
  LDFLAGS += -s
endif

LDFLAGS += `sdl2-config --libs`
CFLAGS += -DSDL2 `sdl2-config --cflags`

all: $(EXECUTABLE)

debug:
	make BUILD=debug

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@
	
.c.o:
	$(CC) $(CFLAGS) $< -o $@

# Add header dependencies here
render.o : render.c sdl/pocadv.h mode7.h obj.h
mode7.o : mode7.c bmp.h mode7.h obj.h
bmp.o : bmp.c bmp.h
pocadv.o : pocadv.c sdl/pocadv.h bmp.h app.h
obj.o : obj.c obj.h

.PHONY : clean

clean:
	-rm -f sdl/*.o *.o $(EXECUTABLE) $(EXECUTABLE).exe
	-rm -f *.log dump.txt

DISTFILE=mode7.zip
dist: clean
	zip -r $(DISTFILE) *
