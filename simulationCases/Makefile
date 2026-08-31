# Optional local environment. If absent, BASILISK must already be set.
ifneq ($(wildcard ../.project_config),)
  BASILISK := $(shell . ../.project_config && echo $$BASILISK)
endif

ifeq ($(BASILISK),)
  $(error BASILISK is unset. Export it or create ../.project_config)
endif

CFLAGS += -O2 -disable-dimensions
CFLAGS += -I$(PWD)/src-local -I$(PWD)/../src-local

include $(BASILISK)/Makefile.defs
