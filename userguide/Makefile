all: doc clean

doc:
	pdflatex djles.tex
	bibtex djles.aux
	pdflatex djles.tex
	pdflatex djles.tex

clean:
	rm -rf *.aux *.bbl *.blg *.log *.out

clobber: clean
	rm -rf *.pdf

