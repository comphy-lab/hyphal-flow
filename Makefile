QCC ?= $(shell command -v qcc 2>/dev/null || { test -n "$$BASILISK" && test -x "$$BASILISK/qcc" && printf '%s\n' "$$BASILISK/qcc"; })
BASILISK ?=
BASILISK_API_HEADER := $(shell qcc_dir="$$(dirname "$(QCC)")"; \
	for header in "$(BASILISK)/grid/multigrid-common.h" \
	              "$(BASILISK)/src/grid/multigrid-common.h" \
	              "$$qcc_dir/grid/multigrid-common.h"; do \
		test -f "$$header" && { printf '%s' "$$header"; break; }; \
	done)
BASILISK_API_FLAG := $(shell if test -z "$(BASILISK_API_HEADER)"; then \
	printf '%s' 'BASILISK_HEADER_NOT_FOUND'; \
	elif ! grep -q 'void set_prolongation' "$(BASILISK_API_HEADER)"; then \
	printf '%s' '-DHYPHAL_LEGACY_BASILISK=1'; \
	fi)
BUILD_DIR ?= build
SOURCE := simulationCases/hyphal-flow.c
TARGET := $(BUILD_DIR)/hyphal-flow
HEADERS := $(wildcard src-local/*.h)

.PHONY: all compile-check help

all: $(TARGET)

$(TARGET): $(SOURCE) $(HEADERS)
	@test -n "$(QCC)" || { echo "qcc not found in PATH" >&2; exit 2; }
	@test "$(BASILISK_API_FLAG)" != "BASILISK_HEADER_NOT_FOUND" || { \
		echo "cannot locate Basilisk grid headers; set BASILISK to the source tree" >&2; \
		exit 2; \
	}
	@mkdir -p "$(BUILD_DIR)"
	"$(QCC)" -I"$(CURDIR)/src-local" -Wall -O2 -disable-dimensions $(BASILISK_API_FLAG) \
		"$(SOURCE)" -o "$@" -lm

compile-check:
	@output_root="$$(mktemp -d "$${TMPDIR:-/tmp}/hyphal-compile-check.XXXXXX")"; \
	trap 'rm -rf "$$output_root"' EXIT; \
	bash runSimulation.sh default.params --compile-only \
		--output-root "$$output_root"

help:
	@echo "make             Compile the canonical simulation into build/"
	@echo "make compile-check  Compile default.params through runSimulation.sh --compile-only"
