(TeX-add-style-hook "tutorial"
 (function
  (lambda ()
    (LaTeX-add-bibitems
     "PAW")
    (TeX-add-symbols
     '("mbax" 1)
     '("vdefault" 1)
     '("vrules" 1)
     '("vformat" 1)
     '("vdescr" 1)
     '("key" 1)
     '("bdescr" 1)
     '("brules" 1)
     '("block" 1))
    (TeX-run-style-hooks
     "shadow"
     "color"
     "graphics"
     "ifthen"
     "times"
     "latex2e"
     "art12"
     "article"
     "final"
     "12pt"))))

