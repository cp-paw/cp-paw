(TeX-add-style-hook "newmanual"
 (function
  (lambda ()
    (LaTeX-add-bibitems
     "Pdcat"
     "PAW"
     "KohnSham"
     "DFTBenchmarks"
     "CP"
     "CP2"
     "Pastore"
     "Madden"
     "Decouple"
     "MPI"
     "ESSL"
     "Constants"
     "Dataexplorer"
     "XMGR"
     "Cerius2"
     "PerdewZunger"
     "CeperleyAlder"
     "GGA91"
     "PW91L"
     "Becke88"
     "Perdew86"
     "PBE"
     "Verlet"
     "Acceleration"
     "Adiabaticity"
     "Nose"
     "Constraints"
     "HinSi"
     "UFF"
     "FreeSoftwareFoundation"
     "LDA+U"
     "MPEGPLAY"
     "MPEGENCODE")
    (LaTeX-add-labels
     "ROOT"
     "fig:h2copdos"
     "constants")
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

