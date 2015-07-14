(TeX-add-style-hook
 "thesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=25mm") ("iwona" "math") ("fontenc" "T1") ("inputenc" "utf8") ("babel" "english") ("biblatex" "backend=biber" "style=alphabetic") ("hyperref" "pdfusetitle" "		pdfauthor={Sudip Sinha}" "		pdfsubject={Masters Thesis}" "		pdfkeywords={derivative, option, pricing, Cox Ross Rubenstein, CRR, Black Scholes, BS}" "		pdfstartview=FitBH" "		pdfpagelayout=OneColumn" "		bookmarks=true" "unicode=true" "colorlinks=true" "linktocpage=false" "			") ("hypcap" "figure" "table")))
   (TeX-run-style-hooks
    "geometry"
    "iwona"
    "eulervm"
    "fontenc"
    "inputenc"
    "babel"
    "csquotes"
    "graphicx"
    "booktabs"
    "parskip"
    "biblatex"
    "pgf"
    "tikz"
    "mathrsfs"
    "color"
    "titlesec"
    "hyperref"
    "hypcap")
   (LaTeX-add-environments
    "thm"
    "crr"
    "prp"
    "lmm"
    "dfn"
    "rem"
    "eg")
   (LaTeX-add-bibliographies)))

