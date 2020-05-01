rnw:
	R -e "library(knitr); knit2pdf('$(FILE)')"

clean:
	\rm -rf *.log *.nav *.aux *.tex *.out *.snm *.toc *.vrb *.bbl ._*

