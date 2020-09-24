LAKES = $(wildcard analysis/data/*_*)
INPUT_TXT = $(patsubst analysis/data/%, analysis/data/%/input.txt, $(LAKES))

test: 
	echo $(INPUT_TXT)
	
all: $(INPUT_TXT)

$(INPUT_TXT): 
	@Rscript -e "ii <- '$(shell dirname $@)'; source('analysis/scripts/01_data_merge.R')"

clean:
	-rm $(INPUT_TXT)