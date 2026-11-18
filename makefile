PROJECT = prime_power_notes
PDFVIEWER ?= pdfterm

.PHONY: default display clean

default:
	latexmk $(PROJECT).tex

display: default
	$(PDFVIEWER) obj/$(PROJECT).pdf

clean:
	latexmk -C -bibtex-cond1 $(PROJECT).tex
