fn=soiliceDocR1

$(fn).pdf: $(fn).tex
	pdflatex $(fn).tex
	pdflatex $(fn).tex
	rm -f $(fn).aux
	rm -f $(fn).log
	rm -f $(fn).out
	rm -f $(fn).synctex.gz
	open $(fn).pdf

.PHONY: clean

clean:
	rm -f $(fn).aux
	rm -f $(fn).log
	rm -f $(fn).out
	rm -f $(fn).synctex.gz
	rm -f $(fn).fdb_latexmk
	rm -f $(fn).fls
	rm -f $(fn).bbl
	rm -f $(fn).blg
