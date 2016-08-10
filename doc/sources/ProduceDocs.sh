#!/bin/bash

for ((index=0; index<2; index++))
do
    pdflatex CC_observables_def.tex
    pdflatex Complex_DGLAP.tex
    pdflatex DIS.tex
    pdflatex Evolution_code.tex
    pdflatex intrinsic_charm.tex
    pdflatex Lagrange_derivative.tex
    pdflatex QCD_QED_common_basis.tex
    pdflatex TruncatedSolution.tex
    pdflatex manual.tex
    pdflatex running_mass.tex
    pdflatex Luminosities.tex
    pdflatex matching_conditions.tex
done

mv *.pdf ../pdfs/.
rm *.aux *.log *.out *.toc *.idx
