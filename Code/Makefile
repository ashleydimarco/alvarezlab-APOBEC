rnw:
	R -e "library(knitr); knit2pdf('$(FILE)')"

rmd:
	Rscript -e "library(knitr); rmarkdown::render('$(FILE)', 'pdf_document')"

rmd_html:
	Rscript -e "library(knitr); rmarkdown::render('$(FILE)', 'html_document')"
	
clean:
	\rm -rf *.log *.nav *.aux *.tex *.out *.snm *.toc *.vrb *.bbl
