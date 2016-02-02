default: ../formulas.png

../formulas.png: formulas.tex
	pdflatex $<
	convert -density 100 -trim formulas.pdf $@


