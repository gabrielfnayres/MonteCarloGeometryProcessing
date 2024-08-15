CC = g++
CFLAGS = -Wall -Wextra
TARGET = montecarlo
SRCS = main.cpp src/monte_carlo.cpp includes/geometry.h
OBJS = $(SRCS:.c=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)