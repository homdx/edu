(TeX-add-style-hook
 "tr_asian_singularpoints"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "hmarginratio=1:1" "top=32mm" "margin=20mm") ("fontenc" "T1") ("babel" "english")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "geometry"
    "lmodern"
    "fontenc"
    "babel"
    "amsmath"
    "amsfonts"
    "amsthm"
    "bm"
    "graphicx"
    "booktabs"
    "hyperref"
    "fancyhdr")
   (TeX-add-symbols
    '("horrule" 1))
   (LaTeX-add-labels
    "eq:arithmeticmean")
   (LaTeX-add-environments
    "thm"
    "corr"
    "lem"
    "prop"
    "rem"
    "defn")))

