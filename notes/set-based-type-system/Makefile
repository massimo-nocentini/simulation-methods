
pdf:
	lualatex -shell-escape article.tex

world:
	lualatex -shell-escape article.tex
	bibtex article
	lualatex -shell-escape article.tex
	lualatex -shell-escape article.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.pdf

jupyter:
	PYTHONPATH=../../src/ jupyter notebook
