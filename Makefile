SHELL := /bin/bash
FOLDERS = $(wildcard analysis/data/*_*)
LAKES = $(shell for lake in $(FOLDERS); do basename $$lake; done)
INPUT_TXT = $(patsubst analysis/data/%, analysis/data/%/input.txt, $(FOLDERS))
# If we want to have ids appended to data files do below:
# TEMPERATURES := $(foreach lake, $(LAKES), analysis/data/$(lake)/$(lake)_temperatures.csv)
# Else
TEMPERATURES := $(foreach lake, $(LAKES), analysis/data/$(lake)/temperatures.csv)
METEOROLOGY := $(foreach lake, $(LAKES), analysis/data/$(lake)/NLDAS.csv)
WATERQUALITY := $(foreach lake, $(LAKES), analysis/data/$(lake)/wq_data.csv)
NML := $(foreach lake, $(LAKES), analysis/data/$(lake)/nhdhr.nml)

test: 
	@echo $(INPUT_TXT)
	
all: $(INPUT_TXT)

%/input.txt: analysis/scripts/01_data_merge.R %/temperatures.csv  %/NLDAS.csv %/wq_data.csv %/nhdhr.nml
	@Rscript -e "ii <- '$(shell dirname $@)'; source('analysis/scripts/01_data_merge.R')"

clean:
	-rm $(INPUT_TXT)
	