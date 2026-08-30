QCC ?= $(shell command -v qcc 2>/dev/null || { test -n "$$BASILISK" && test -x "$$BASILISK/qcc" && printf '%s\n' "$$BASILISK/qcc"; })
BUILD_DIR ?= build
SOURCE := simulationCases/hyphal-flow.c
TARGET := $(BUILD_DIR)/hyphal-flow
HEADERS := $(wildcard src-local/*.h)

.PHONY: all compile-smoke help

all: $(TARGET)

$(TARGET): $(SOURCE) $(HEADERS)
	@test -n "$(QCC)" || { echo "qcc not found in PATH" >&2; exit 2; }
	@mkdir -p "$(BUILD_DIR)"
	"$(QCC)" -I"$(CURDIR)/src-local" -Wall -O2 -disable-dimensions \
		"$(SOURCE)" -o "$@" -lm

compile-smoke:
	@output_root="$$(mktemp -d "$${TMPDIR:-/tmp}/hyphal-compile-smoke.XXXXXX")"; \
	trap 'rm -rf "$$output_root"' EXIT; \
	bash runSimulation.sh smoke/full.params --compile-only \
		--output-root "$$output_root"

help:
	@echo "make             Compile the canonical simulation into build/"
	@echo "make compile-smoke  Exercise the clean-clone runner compile path"
