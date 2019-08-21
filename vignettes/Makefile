all: pdf 

pdf: DClusterm.pdf

%.pdf: %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*
	pdflatex $*
	#rm DClusterm-0*pdf

clean:
	rm *tex *aux *log *bbl *blg *out DClusterm-*.pdf jss.bst jss.cls
	
