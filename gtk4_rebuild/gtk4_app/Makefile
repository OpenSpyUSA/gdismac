APP := gdis-gtk4
BUILD_DIR := build
SRC_DIR := src
TOOL_DIR := tools
MACOS_DIR := macos
LEGACY_SGINFO_DIR := ../legacy_snapshot/src
APP_BUNDLE_NAME := GDIS.app
APP_BUNDLE := $(BUILD_DIR)/$(APP_BUNDLE_NAME)
APP_BUNDLE_CONTENTS := $(APP_BUNDLE)/Contents
APP_BUNDLE_MACOS := $(APP_BUNDLE_CONTENTS)/MacOS
APP_BUNDLE_RESOURCES := $(APP_BUNDLE_CONTENTS)/Resources
APP_BUNDLE_PLIST_SRC := $(MACOS_DIR)/Info.plist
APP_BUNDLE_LAUNCHER := $(APP_BUNDLE_MACOS)/$(APP)
APP_BUNDLE_BINARY := $(APP_BUNDLE_MACOS)/$(APP)-bin
DIST_DIR := ../../dist
DIST_BUNDLE_NAME := GDIS.app
DIST_BUNDLE := $(DIST_DIR)/$(DIST_BUNDLE_NAME)
PORTABLE_RELEASE_BUILDER := $(TOOL_DIR)/build_macos_portable_release.py
SOURCES := \
	$(SRC_DIR)/main.c \
	$(SRC_DIR)/gdis_elements.c \
	$(SRC_DIR)/gdis_gtk4_window.c \
	$(SRC_DIR)/gdis_isosurface.c \
	$(SRC_DIR)/gdis_model.c \
	$(SRC_DIR)/gdis_sgcore.c \
	$(SRC_DIR)/gdis_restoration.c \
	$(SRC_DIR)/gdis_reports.c \
	$(SRC_DIR)/gdis_legacy_map.c
LEGACY_SOURCES := \
	$(LEGACY_SGINFO_DIR)/sgclib.c \
	$(LEGACY_SGINFO_DIR)/sgio.c
TOOL_SOURCES := \
	$(TOOL_DIR)/gdis_model_probe.c \
	$(TOOL_DIR)/gdis_model_bond_audit.c

GTK4_CFLAGS := $(shell pkg-config --cflags gtk4 2>/dev/null)
GTK4_LIBS := $(shell pkg-config --libs gtk4 2>/dev/null)

CFLAGS += -std=c11 -Wall -Wextra -Wpedantic -Wno-extra-semi -g -MMD -MP -I$(LEGACY_SGINFO_DIR)
OBJCFLAGS += -Wall -Wextra -g -MMD -MP
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
OBJECTS += $(patsubst $(LEGACY_SGINFO_DIR)/%.c,$(BUILD_DIR)/legacy_%.o,$(LEGACY_SOURCES))
TOOL_OBJECTS := $(patsubst $(TOOL_DIR)/%.c,$(BUILD_DIR)/%.tool.o,$(TOOL_SOURCES))
DEPS := $(OBJECTS:.o=.d) $(TOOL_OBJECTS:.o=.d)
TOOL_LINK_OBJECTS := \
	$(BUILD_DIR)/gdis_model.o \
	$(BUILD_DIR)/gdis_elements.o \
	$(BUILD_DIR)/gdis_sgcore.o \
	$(BUILD_DIR)/legacy_sgclib.o \
	$(BUILD_DIR)/legacy_sgio.o
MODEL_PROBE := $(BUILD_DIR)/gdis-model-probe
BOND_AUDIT := $(BUILD_DIR)/gdis-model-bond-audit

.PHONY: all clean run doctor tools app open-app dist-app open-dist-app portable-release export-png

-include $(DEPS)

all: doctor $(BUILD_DIR)/$(APP) tools

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

tools: $(MODEL_PROBE) $(BOND_AUDIT)

$(MODEL_PROBE): $(BUILD_DIR)/gdis_model_probe.tool.o $(TOOL_LINK_OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BOND_AUDIT): $(BUILD_DIR)/gdis_model_bond_audit.tool.o $(TOOL_LINK_OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(GTK4_CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.m
	@mkdir -p $(BUILD_DIR)
	$(CC) $(OBJCFLAGS) $(GTK4_CFLAGS) -c $< -o $@

$(BUILD_DIR)/legacy_%.o: $(LEGACY_SGINFO_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(GTK4_CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.tool.o: $(TOOL_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(GTK4_CFLAGS) -c $< -o $@

run: all
	./$(BUILD_DIR)/$(APP) $(ARGS)

export-png: all
	@: $${MODEL:?'Set MODEL=/path/to/model'}
	@: $${OUT:?'Set OUT=/path/to/output.png'}
	GDIS_GTK4_NON_UNIQUE=1 \
	GDIS_GTK4_FORCE_STARTUP_PATHS=1 \
	GDIS_GTK4_STARTUP_MODEL_PATH="$(MODEL)" \
	GDIS_GTK4_STARTUP_EXPORT_PNG="$(OUT)" \
	GDIS_GTK4_STARTUP_EXPORT_SIZE="$(if $(SIZE),$(SIZE),1360x852)" \
	GDIS_GTK4_STARTUP_CAPTURE_ISO="$(if $(ISO),$(ISO),0)" \
	GDIS_GTK4_STARTUP_ISO_MODE="$(MODE)" \
	GDIS_GTK4_STARTUP_ISO_COLOR="$(COLOR)" \
	GDIS_GTK4_STARTUP_EXPORT_AND_QUIT=1 \
	./$(BUILD_DIR)/$(APP) "$(MODEL)"

app: all
ifeq ($(UNAME_S),Darwin)
	@mkdir -p "$(APP_BUNDLE_MACOS)" "$(APP_BUNDLE_RESOURCES)"
	@cp "$(APP_BUNDLE_PLIST_SRC)" "$(APP_BUNDLE_CONTENTS)/Info.plist"
	@cp "$(MACOS_DIR)/gdis-gtk4-launcher.sh" "$(APP_BUNDLE_LAUNCHER)"
	@chmod +x "$(APP_BUNDLE_LAUNCHER)"
	@cp "$(BUILD_DIR)/$(APP)" "$(APP_BUNDLE_BINARY)"
	@chmod +x "$(APP_BUNDLE_BINARY)"
	@if command -v codesign >/dev/null 2>&1; then \
		codesign --force --deep --sign - "$(APP_BUNDLE)" >/dev/null; \
	fi
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
	@if command -v codesign >/dev/null 2>&1; then \
		codesign --force --deep --sign - "$(DIST_BUNDLE)" >/dev/null; \
	fi
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

portable-release: all
ifeq ($(UNAME_S),Darwin)
	python3 "$(PORTABLE_RELEASE_BUILDER)" --skip-build
else
	@echo "portable-release is only available on macOS."
endif

clean:
	rm -rf $(BUILD_DIR)
