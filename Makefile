################################################################################
# This is a general-purpose Makefile for building C/C++ projects.
#
# Some features contained within:
#   * automatic dependency tracking, header files are parsed out of source
#     and targets generated so that when they change, only the needed binaries
#     are re-compiled.
#   * distribution tarball generation (make dist)
#   * header file verification (make hcheck)
################################################################################

# force the use of bash to avoid failing on systems like Ubuntu
SHELL := /bin/bash

# helper function for printing information.  note that filter a,b is Make-ese for
# a == b in this context.  So if $1 == $2, print $3.
define printif
	$(if $(filter $1,$2), $(info $3))
endef


################################################################################
# configure compilation environment
################################################################################

# look for gcc >= 4.8.5
ICC = $(shell command -v icc 2> /dev/null)
ifdef ICC
		ICC_VER  = $(shell echo `icc -dumpversion | awk -F. '{ printf "%02i%02i%02i\n", $$1, $$2, $$3}'`)
		ICC_OK = $(shell expr $(ICC_VER) \>= 130000)
else
		ICC_OK = 0
endif
GCC_OK = $(shell expr `gcc -dumpversion | awk -F. '{ printf "%02i%02i%02i\n", $$1, $$2, $$3}'` \>= 040805)
ifeq "$(ICC_OK)" "1"
    CXX := icc
    CC  := icc
else ifeq "$(GCC_OK)" "1"
    CXX := g++
    CC  := gcc
else
    $(error gcc >= 4.8.5 or intel >= 13 is required to compile)
endif

MKLPATH  = $(MKLROOT)/lib/intel64
# build link list for MKL libraries, note gcc requires a different threading library
MKLLDOPT = -Wl,--start-group $(MKLPATH)/libmkl_intel_ilp64.a $(MKLPATH)/libmkl_core.a
ifeq ($(CXX),g++)
		MKLLDOPT += $(MKLPATH)/libmkl_gnu_thread.a
else
		MKLLDOPT += $(MKLPATH)/libmkl_intel_thread.a
endif
MKLLDOPT += -Wl,--end-group

# compiler configuration
override CXXOPTS := -std=c++11 -O3 -g3 -Wall -Wextra -fno-omit-frame-pointer -fopenmp -pthread $(CXXOPTS)
override LDOPTS  := -Wl,--as-needed $(MKLLDOPT) $(LDOPTS) -lrt -Wl,--no-as-needed -ldl

# add coverage flags if requested
ifeq ($(COVERAGE),1)
    override CXXOPTS += -fprofile-arcs -ftest-coverage
endif

# add debug flags if requested
ifeq ($(DEBUG),1)
    override CXXOPTS += -O0
endif

# add profiling flags if requested
ifeq ($(PROFILE),1)
    override CXXOPTS += -pg
endif

# supress warning #2196 on intel as it's nonsense, and inlining warnings
ifeq ($(CXX),icc)
    override CXXOPTS += -wd2196 -diag-disable remark
endif

# enable fast math on GCC by default (note this is also default behavior for Intel), also disable
# stack protection as we don't need it and it has a small performance cost
ifeq ($(CXX),g++)
    override CXXOPTS += -ffast-math
    override CXXOPTS += -fno-stack-protector
endif

# determine architecture
ARCH = $(shell uname -m)

# export variables for sub-shells
export CC CXX ARCH VERBOSE

# print configuration information
$(info configuration)
$(call printif,1,             1,  · $(CXX) $(shell $(CXX) -dumpversion))
$(call printif,1,             1,  · $(ARCH))
$(call printif,$(DEBUG),      1,  · debug)
$(call printif,$(COVERAGE),   1,  · enabling coverage support)
$(call printif,$(PROFILE),    1,  · enabling profiling support)
$(info ${EMPTY})

# silent by default
ifndef VERBOSE
.SILENT:
endif

################################################################################
# configure project
################################################################################

# define project name and help.  This will tell us the tarball to generate
# and provide help text to print when "make help" is run
PROJECT_NAME := klt
PROJECT_HELP := "KLT library"

# grab binaries to build
PROGRAMS := $(patsubst src/%.cc, bin/%, $(shell echo `ls src/*.cc`))

# all sources for docsnip
SOURCES  := $(shell find src/ test/ -iname "*.h" -o -iname "*.cc" -o -iname "*.py" | grep -v '\#')


################################################################################
# Define targets
################################################################################
.PHONY: all
all: $(PROGRAMS) ## (default target) build everything

# print help string
COLWIDTH = 15

.PHONY: help
help:                          ## print this message and exit
	@echo $(PROJECT_HELP)
	@echo
	@echo "Useful command-line variables:"
	@echo "  VERBOSE=1  -- enable verbose mode, print compilation commands"
	@echo "  COVERAGE=1 -- add flags for coverage information"
	@echo "  DEBUG=1    -- enable debug build (no optimization)"
	@echo "  PROFILE=1  -- enable flags for profiling"
	@echo
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:[[:space:]]*.*?## / {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST) | sort


# build and run test harness
.PHONY: test
test: bin/test               ## build and run tests
	@echo " ▸ running doctest tests"
	@echo
	@bin/test


# General binary target for c++
bin/% : src/%.cc
	@rm -f $@
	@$(CXX) $(CXXOPTS) $<
	@echo " ▸   [BIN] $(basename $(notdir $@))"
	$(CXX) $(CXXOPTS) $< -o $@ $(LDOPTS)


# source distribution target
RELNAME = $(shell if [ `git describe --tags >& /dev/null` ]; then git describe --tags; else git rev-parse --short HEAD; fi)
.PHONY: srcdist
srcdist: $(PROJECT_NAME)-$(RELNAME)-src.tar.xz   ## build source distribution tarball

$(PROJECT_NAME)-$(RELNAME)-src.tar.xz:
	@echo " ▸ building $@"
	@rm -f $@
	@echo "   ▸ cloning to /tmp"
	@rm -rf /tmp/$(PROJECT_NAME)
	@git clone ./ /tmp/$(PROJECT_NAME) >& /dev/null
	@echo
	@cd /tmp/$(PROJECT_NAME)              && \
	  rm -rf `find . -iname ".git*"`; rm -rf .objs; cd .. && \
	  tar -cf - $(PROJECT_NAME) | xz -9 -c - > $@
	@mv /tmp/$@ ./
	@rm -rf /tmp/$(PROJECT_NAME)


# binary distribution target
.PHONY: dist
dist: $(PROJECT_NAME)-$(RELNAME).tar.xz ## build binary distribution tarball

$(PROJECT_NAME)-$(RELNAME).tar.xz:
	@echo " ▸ building $@"
	@rm -f $@
	@echo "   ▸ cloning to /tmp"
	@rm -rf /tmp/$(PROJECT_NAME)
	@git clone ./ /tmp/$(PROJECT_NAME) >& /dev/null
#@echo "   ▸ building static binaries"
	@echo
	@cd /tmp/$(PROJECT_NAME)              && \
	  $(MAKE) --no-print-directory && \
	  rm -rf `find . -iname ".git*"`; rm -rf .objs; cd .. && \
	  tar -cf - $(PROJECT_NAME) | xz -9 -c - > $@
	@mv /tmp/$@ ./
	@rm -rf /tmp/$(PROJECT_NAME)


# cleanup targets
.PHONY: clean
clean:                         ## clean normal build detritus
	@rm -f $(PROGRAMS)
	@rm -f bin/test

.PHONY: remake
remake: clean all              ## clean everything and rebuild


################################################################################
# check header independence
################################################################################

# hcheck attempts to compile each header file as a standalone binary
# this is useful for ensuring consistency, and that each headerfile properly
# includes all its required dependencies so it can easily be sliced out.
HEADERSRC = $(shell echo `find src/ -iname "*.h" | sort`)
HEADEROBJ = $(patsubst %.h, %.h.o, $(HEADERSRC))

%.h.cc: %.h
	@grep -v "#pragma once" $< > $@

%.h.o: %.h.cc
	@echo " ▸ [HCHECK] $(basename $*.h)"
	$(CXX) $(CXXOPTS) -c $< -o $@
	@rm $@ $<

.PHONY: hcheck
hcheck: $(HEADEROBJ)           ## run header check to verify consistency and isolation of headers


################################################################################
# boilerplate
################################################################################

# include any local make modifications (for testing normally)
-include make.local
