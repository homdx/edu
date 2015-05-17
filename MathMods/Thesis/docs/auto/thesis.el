(TeX-add-style-hook
 "thesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=25mm") ("iwona" "math") ("fontenc" "T1") ("inputenc" "utf8") ("babel" "english") ("microtype" "protrusion=true" "expansion=true" "final" "babel") ("hyperref" "pdfusetitle" "pdfauthor={Sudip Sinha}" "pdfsubject={Masters Thesis}" "pdfkeywords={derivative, option, pricing, Cox Ross Rubenstein, CRR, Black Scholes, BS}" "pdfstartview=FitBH" "pdfpagelayout=OneColumn" "bookmarks=true" "unicode=true" "colorlinks=false" "") ("hypcap" "figure" "table")))
   (TeX-run-style-hooks
    "geometry"
    "iwona"
    "eulervm"
    "fontenc"
    "inputenc"
    "babel"
    "amsmath"
    "amsfonts"
    "amssymb"
    "graphicx"
    "microtype"
    "color"
    "titlesec"
    "lastpage"
    "fancyhdr"
    "hyperref"
    "hypcap"
    "booktabs")
   (TeX-add-symbols
    "LaTeXtitle")
   (LaTeX-add-environments
    "thm"
    "crr"
    "prp"
    "lmm"
    "dfn"
    "rem")))

