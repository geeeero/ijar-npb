# how to make eps files for submission of figures if necessary

# first run this
library(grDevices)
# then use the (commented out) lines with "cairo_ps(...)" instead of the "pdf(...)" lines
# (so far only made for figures in file "pdc-illu.R")
# NB: simply using postscript(...) does not work for transparency

# figures "singleprior-pdc2" (Fig 1), "betaset-binomset1" (Fig 3) and "paramsets" (Fig 2)
# are in file "pdc-illu.R".

# Figure 4 is pure tikz, so no pdf/eps file

# figures "bridge-fittingfailures" (Fig 5), "bridge-earlyfailures" (Fig 6)
# and "bridge-latefailures" (Fig 7)
# are in file "bridgesystem.R".

# figure "brakingsystem-2" is in file "brakingsystem.R"

# figures "need-to-opt-x", x = 1:8, are in file "need-to-opt.R"

#