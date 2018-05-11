
dmcDoc.pdf: dmcDoc.tex dmcDoc.bib
	latexmk -f -pdf -pdflatex="pdflatex -interaction=nonstopmode" -bibtex -use-make dmcDoc.tex

clean:
	latexmk -C

