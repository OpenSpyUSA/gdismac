APP := gdis-gtk4
BUILD_DIR := build
SRC_DIR := src
SOURCES := \
	$(SRC_DIR)/main.c \
	$(SRC_DIR)/gdis_elements.c \
	$(SRC_DIR)/gdis_gtk4_window.c \
	$(SRC_DIR)/gdis_isosurface.c \
	$(SRC_DIR)/gdis_model.c \
	$(SRC_DIR)/gdis_restoration.c \
	$(SRC_DIR)/gdis_reports.c \
	$(SRC_DIR)/gdis_legacy_map.c

GTK4_CFLAGS := $(shell pkg-config --cflags gtk4 2>/dev/null)
GTK4_LIBS := $(shell pkg-config --libs gtk4 2>/dev/null)

CFLAGS += -std=c11 -Wall -Wextra -Wpedantic -Wno-extra-semi -g
OBJCFLAGS += -Wall -Wextra -g
LDFLAGS +=
LDLIBS += $(GTK4_LIBS)
ARGS ?=

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
SOURCES += $(SRC_DIR)/gdis_macos_menu.m
LDLIBS += -framework AppKit -framework Foundation
endif

OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(filter %.c,$(SOURCES)))
OBJECTS += $(patsubst $(SRC_DIR)/%.m,$(BUILD_DIR)/%.o,$(filter %.m,$(SOURCES)))

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

$(BUILD_DIR)/$(APP): $(OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(GTK4_CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.m
	@mkdir -p $(BUILD_DIR)
	$(CC) $(OBJCFLAGS) $(GTK4_CFLAGS) -c $< -o $@

run: all
	./$(BUILD_DIR)/$(APP) $(ARGS)

clean:
	rm -rf $(BUILD_DIR)
