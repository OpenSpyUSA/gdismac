APP := gdis-gtk4
BUILD_DIR := build
SRC_DIR := src
SOURCES := \
	$(SRC_DIR)/main.c \
	$(SRC_DIR)/gdis_gtk4_window.c \
	$(SRC_DIR)/gdis_legacy_map.c

GTK4_CFLAGS := $(shell pkg-config --cflags gtk4 2>/dev/null)
GTK4_LIBS := $(shell pkg-config --libs gtk4 2>/dev/null)

CFLAGS += -std=c11 -Wall -Wextra -Wpedantic -g

.PHONY: all clean run doctor

all: doctor $(BUILD_DIR)/$(APP)

doctor:
	@if ! command -v pkg-config >/dev/null 2>&1; then \
		echo "pkg-config is required."; \
		echo "Install it first, then install GTK4."; \
		exit 1; \
	fi
	@if ! pkg-config --exists gtk4; then \
		echo "GTK4 development files were not found."; \
		echo "On macOS with Homebrew: brew install gtk4 pkgconf"; \
		exit 1; \
	fi

$(BUILD_DIR)/$(APP): $(SOURCES)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(GTK4_CFLAGS) $(SOURCES) -o $@ $(GTK4_LIBS)

run: all
	./$(BUILD_DIR)/$(APP)

clean:
	rm -rf $(BUILD_DIR)
