TEX = mlegp.tex gp.tex gp_ex.tex SA.tex SVD.tex overview.tex
RNW = $(TEX:.tex=.Rnw)
PDF = $(TEX:.tex=.pdf)
PS = $(TEX:.tex=.ps)

mlegp.dvi: ${TEX}
	latex mlegp.tex 
	latex mlegp.tex 
	latex mlegp.tex 

### update tex files ###
mlegp.tex: mlegp.Rnw
	R CMD Sweave mlegp.Rnw

gp.tex: gp.Rnw
	R CMD Sweave gp.Rnw

gp_ex.tex: gp_ex.Rnw
	R CMD Sweave gp_ex.Rnw

SA.tex: SA.Rnw
	R CMD Sweave SA.Rnw

SVD.tex: SVD.Rnw
	R CMD Sweave SVD.Rnw

overview.tex: overview.Rnw
	R CMD Sweave overview.Rnw

### update ps files ##
#mlegp.ps: mlegp.dvi
#	dvips -Pwww -t letter mlegp.dvi -o mlegp.ps



clean:
	rm -f *.aux *.log *.bbl *.blg *.eps ${TEX} *.ps *.dvi 



