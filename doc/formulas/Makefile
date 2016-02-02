default: ../formulas.png

../formulas.png: formulas.tex
	pdflatex $<
	convert -density 150 -trim formulas.pdf $@


