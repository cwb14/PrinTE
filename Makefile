# PrinTE build and check targets. Nothing here needs root.

PREFIX  ?= $(CURDIR)
BIN     := $(PREFIX)/bin
MUTATOR := $(BIN)/ltr_mutator
SRC     := src/printe/cpp/ltr_mutator.cpp
VERSION := $(shell sed -n 's/^__version__ = "\(.*\)"/\1/p' src/printe/__init__.py)
CACHE   := $(if $(XDG_CACHE_HOME),$(XDG_CACHE_HOME),$(HOME)/.cache)/printe

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  # Apple clang has no -fopenmp driver flag; it needs libomp from Homebrew and a
  # static link against it.
  OMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
  CXX        ?= clang++
  OMPFLAGS   := -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
  OMPLIBS    := $(OMP_PREFIX)/lib/libomp.a
else
  CXX      ?= g++
  OMPFLAGS := -fopenmp
  OMPLIBS  :=
endif

.PHONY: all ltr-mutator test test-fast lint fetch-data docker apptainer clean help

help:
	@echo "PrinTE $(VERSION)"
	@echo "  make ltr-mutator   build the point-mutation binary into $(BIN)"
	@echo "  make test          run the full test suite"
	@echo "  make test-fast     skip the end-to-end tests"
	@echo "  make lint          ruff + shellcheck"
	@echo "  make fetch-data    download the optional LTR-RT exemplar database"
	@echo "  make docker        build the container image"
	@echo "  make apptainer     build a .sif from the published image"
	@echo "  make clean         remove build output and caches"

all: ltr-mutator

ltr-mutator: $(MUTATOR)

$(MUTATOR): $(SRC)
	@mkdir -p $(BIN)
ifeq ($(UNAME),Darwin)
	@test -n "$(OMP_PREFIX)" || { echo "libomp not found. Run: brew install libomp"; exit 1; }
endif
	$(CXX) -std=c++17 -O3 $(OMPFLAGS) $< -o $@ $(OMPLIBS)

test: ltr-mutator
	PYTHONPATH=src python -m pytest tests -v

test-fast: ltr-mutator
	PYTHONPATH=src python -m pytest tests -v -m "not slow"

lint:
	ruff check src tests
	shellcheck -S warning PrinTE.sh scripts/*.sh

fetch-data:
	@mkdir -p $(CACHE)/data
	curl -fsSL -o $(CACHE)/data/ltr-db.fa.gz \
	  https://github.com/cwb14/PrinTE/releases/download/v$(VERSION)/ltr-db.fa.gz
	@echo "Fetched into $(CACHE)/data"

docker:
	docker build -t printe:$(VERSION) -t printe:latest .

apptainer:
	apptainer build printe_$(VERSION).sif containers/printe.def

clean:
	rm -rf $(BIN) build dist src/*.egg-info .pytest_cache .ruff_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
