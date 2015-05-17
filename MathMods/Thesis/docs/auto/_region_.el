(TeX-add-style-hook
 "_region_"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "12pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8x") ("babel" "english") ("microtype" "final" "babel") ("hyperref" "backref='section'" "pdftitle='Thesis'" "pdfauthor='SudipSinha'" "pdfsubject='MathematicalFinance'" "pdfkeywords='derivativepricing'" "pdfstartview=FitBH" "pdfpagelayout=OneColumn" "")))
   (TeX-run-style-hooks
    "latex2e"
    "report"
    "rep12"
    "lmodern"
    "inputenc"
    "babel"
    "ucs"
    "amsmath"
    "amsfonts"
    "amssymb"
    "graphicx"
    "microtype"
    "hyperref"
    "booktabs"
    "fancyhdr")))

