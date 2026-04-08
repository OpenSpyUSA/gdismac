APP := gdis-gtk4
BUILD_DIR := build
SRC_DIR := src
MACOS_DIR := macos
APP_BUNDLE_NAME := GDIS.app
APP_BUNDLE := $(BUILD_DIR)/$(APP_BUNDLE_NAME)
APP_BUNDLE_CONTENTS := $(APP_BUNDLE)/Contents
APP_BUNDLE_MACOS := $(APP_BUNDLE_CONTENTS)/MacOS
APP_BUNDLE_RESOURCES := $(APP_BUNDLE_CONTENTS)/Resources
APP_BUNDLE_PLIST_SRC := $(MACOS_DIR)/Info.plist
APP_BUNDLE_LAUNCHER := $(APP_BUNDLE_MACOS)/$(APP)
DIST_DIR := ../../dist
DIST_BUNDLE_NAME := GDIS.app
DIST_BUNDLE := $(DIST_DIR)/$(DIST_BUNDLE_NAME)
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

.PHONY: all clean run doctor app open-app dist-app open-dist-app

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

app: all
ifeq ($(UNAME_S),Darwin)
	@mkdir -p "$(APP_BUNDLE_MACOS)" "$(APP_BUNDLE_RESOURCES)"
	@cp "$(APP_BUNDLE_PLIST_SRC)" "$(APP_BUNDLE_CONTENTS)/Info.plist"
	@cp "$(BUILD_DIR)/$(APP)" "$(APP_BUNDLE_LAUNCHER)"
	@chmod +x "$(APP_BUNDLE_LAUNCHER)"
	@rm -f "$(APP_BUNDLE_MACOS)/$(APP)-bin"
	@echo "Created $(APP_BUNDLE)"
else
	@echo "The app bundle target is only available on macOS."
endif

open-app: app
ifeq ($(UNAME_S),Darwin)
	open "$(APP_BUNDLE)"
else
	@echo "open-app is only available on macOS."
endif

dist-app: app
ifeq ($(UNAME_S),Darwin)
	@mkdir -p "$(DIST_DIR)"
	@rm -rf "$(DIST_BUNDLE)"
	@rm -rf "$(DIST_DIR)/GDIS GTK4.app"
	@cp -R "$(APP_BUNDLE)" "$(DIST_BUNDLE)"
	@echo "Created $(DIST_BUNDLE)"
else
	@echo "dist-app is only available on macOS."
endif

open-dist-app: dist-app
ifeq ($(UNAME_S),Darwin)
	open "$(DIST_BUNDLE)"
else
	@echo "open-dist-app is only available on macOS."
endif

clean:
	rm -rf $(BUILD_DIR)
