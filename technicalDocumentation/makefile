soiliceDoc.pdf: soiliceDoc.tex
	pdflatex soiliceDoc.tex
	pdflatex soiliceDoc.tex
	rm -f soiliceDoc.aux
	rm -r soiliceDoc.log
	rm -f soiliceDoc.out
	rm -f soiliceDoc.synctex.gz
	open soiliceDoc.pdf

.PHONY: clean

clean:
	rm -f soiliceDoc.aux
	rm -f soiliceDoc.log
	rm -f soiliceDoc.out
	rm -f soiliceDoc.synctex.gz
	rm -f soiliceDoc.fdb_latexmk
	rm -f soiliceDoc.fls
	rm -f soiliceDoc.bbl
	rm -f soiliceDoc.blg
