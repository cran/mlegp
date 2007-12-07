TEX = mlegp.tex gp.tex gp_ex.tex SA.tex SVD.tex overview.tex
RNW = $(TEX:.tex=.Rnw)
PDF = $(TEX:.tex=.pdf)
PS = $(TEX:.tex=.ps)

all: ${PDF} clean

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
mlegp.ps: mlegp.dvi
	dvips -Pwww -t letter mlegp.dvi -o mlegp.ps

gp.ps: mlegp.dvi
	dvips -Pwww -pp1-3 -t letter mlegp.dvi -o gp.ps

gp_ex.ps: mlegp.dvi
	dvips -Pwww -pp3-6 -t letter mlegp.dvi -o gp_ex.ps

SA.ps: mlegp.dvi
	dvips -Pwww -pp6-9 -t letter mlegp.dvi -o SA.ps

SVD.ps: mlegp.dvi
	dvips -Pwww -pp10-15 -t letter mlegp.dvi -o SVD.ps

overview.ps: mlegp.dvi
	dvips -Pwww -pp1 -t letter mlegp.dvi -o overview.ps

### update pdf files ###
mlegp.pdf: mlegp.ps
	ps2pdf mlegp.ps

gp.pdf: gp.ps
	ps2pdf gp.ps

gp_ex.pdf: gp_ex.ps
	ps2pdf gp_ex.ps

SA.pdf: SA.ps
	ps2pdf SA.ps

SVD.pdf: SVD.ps
	ps2pdf SVD.ps

overview.pdf: overview.ps
	ps2pdf overview.ps

clean:
	rm -f *.aux *.log *.bbl *.blg *.eps ${TEX} *.ps *.dvi 



