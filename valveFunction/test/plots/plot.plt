#!/bin/gnuplot

set term post enhanced eps color
set output "trans.eps"

alphaOpen = 0.0
beta = pi/10.0

h = 42/1000.0

xf(x) = (x < (alphaOpen + beta)) ? (h/2.0)*sin((pi/2.0) - x) : 1/0 
yf(x) = (x < (alphaOpen + beta)) ? (h/2.0)*sin((pi/2.0) - x) : 1/0 

set table ("lines.dat")
set xrange [0:pi/2.0]
plot xf(x),\
	 yf(x)
unset table


!open trans.eps
