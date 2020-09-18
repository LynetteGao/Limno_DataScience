
all: data

data: analysis/data/143249470_Mendota/input.txt

analysis/data/143249470_Mendota/input.txt: 
	Rscript -e "ii <- '$(shell dirname $@)'; print(ii); source('analysis/scripts/01_data_merge.R')"

clean:
	-rm analysis/data/143249470_Mendota/input.txt
	