CC = g++
CFLAGS = -O2 -Wall
LIBS = -llapacke -lblas

SRC = advection.cpp
TARGET = advection

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

.PHONY: clean

clean:
	rm -f $(TARGET)
