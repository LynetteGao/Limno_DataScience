SHELL := /bin/bash
FOLDERS = $(wildcard analysis/data/*_*)
LAKES = $(shell for lake in $(FOLDERS); do basename $$lake; done)
INPUT_TXT = $(patsubst analysis/data/%, analysis/data/%/input.txt, $(LAKES))
# TEMPERATURES := $(foreach lake, $(LAKES), analysis/data/$(lake)/$(lake)_temperatures.csv)
TEMPERATURES := $(foreach lake, $(LAKES), analysis/data/$(lake)/temperatures.csv)

METEOROLOGY = $(wildcard analysis/data/*/NLDAS*.csv)
WATERQUALITY_ALL = $(wildcard analysis/data/*/wq_data*.csv)
WATERQUALITY_ORIGINAL = $(wildcard analysis/data/*/wq_data*original.csv)
WATERQUALITY = $(filter-out $(WATERQUALITY_ORIGINAL),$(WATERQUALITY_ALL))
NML = $(wildcard analysis/data/*/*.nml)

test: 
	@echo $(TEMPERATURES)
	
all: $(INPUT_TXT)

$(INPUT_TXT): analysis/scripts/01_data_merge.R 
	@Rscript -e "ii <- '$(shell dirname $@)'; source('analysis/scripts/01_data_merge.R')"

clean:
	-rm $(INPUT_TXT)
	