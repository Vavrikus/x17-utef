#!/bin/bash

# Remove LaTeX auxiliary files
rm -f *.aux *.bbl *.blg *.log *.lof *.lot *.fls *.fdb_latexmk \
      *.toc *.synctex.gz *.out *.nav *.snm *.vrb *.bcf *.run.xml

echo "Cleaned LaTeX auxiliary files."
