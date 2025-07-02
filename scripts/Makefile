CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99
LIBS = -lz

# Target executable
TARGET = 1_preprocess_and_trim

# Source files
SOURCES = 1_preprocess_and_trim.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Default target
all: $(TARGET)

# Build the main target
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# Build object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Install target (optional)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Uninstall target (optional)
uninstall:
	rm -f /usr/local/bin/$(TARGET)

# Test target
test: $(TARGET)
	@echo "Testing $(TARGET)..."
	./$(TARGET) --help

# Debug build
debug: CFLAGS += -g -DDEBUG
debug: $(TARGET)

# Profile build
profile: CFLAGS += -pg
profile: $(TARGET)

.PHONY: all clean install uninstall test debug profile 