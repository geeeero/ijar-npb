# to embed fonts in figure pdfs
# (trying to avoid redoing all graphs as eps)
# http://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html

#install.packages("extrafont")
library(extrafont)
font_import()
loadfonts()

embed_fonts("need-to-opt-1.pdf") # does not work

#