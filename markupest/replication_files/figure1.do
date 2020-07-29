f/* Replica of Figure 1 - Markup estimation using Stata: Micro and Macro approaches with markupest by G Rovigatti */
use "${datadir}/mkup_example", clear
/* xtset the data - firms and years */
xtset f_id year
/* set the seed and run the estimation: value added with ACF translog - print results on-screen and correct for first-stage residuals */
set seed 12356
bys nnace: markupest mkupACF_translog, method(dlw) output(ln_va) inputvar(ln_l) free(ln_l) state(ln_k) proxy(ln_m) valueadded prodestopt("poly(3) acf trans") verbose corr
/* plot the graph */
tw (kdensity mkupACF_translog if nnace == 1, lw(medthick) lp(_) lc(ebblue)) (kdensity mkupACF_translog if nnace == 2, lw(medthick) lp(-) lc(maroon)) /*
	*/ (kdensity mkupACF_translog if nnace == 3, lw(medthick) lp(.-.) lc(forest_green)) (kdensity mkupACF_translog if nnace == 4, lw(medthick) lp(-.-) lc(sand)) /*
	*/ (kdensity mkupACF_translog if nnace == 5, lw(medthick) lp(l) lc(navy)) (kdensity mkupACF_translog if nnace == 6, lw(medthick) lp(dot) lc(purple)) /*
	*/ (kdensity mkupACF_translog if nnace == 7, lw(medthick) lp(_) lc(olive_teal)) (kdensity mkupACF_translog if nnace == 8, lw(medthick) lp(-) lc(cyan)) /*
	*/ (kdensity mkupACF_translog if nnace == 9, lw(medthick) lp(l) lc(ltblue)) (kdensity mkupACF_translog if nnace == 10, lw(medthick) lp(dot) lc(mint)) /*
	*/ (kdensity mkupACF_translog if nnace == 11, lw(medthick) lp(_) lc(erose))  /*
	*/ if mkupACF_translog > 0 & mkupACF_translog < 3, ytitle("Density") xtitle("Markup") legend(order( 1 "25" 2 "41" 3 "43" 4 "45" 5 "49" 6 "55" 7 "56" 8 "62" 9 "63" 10 "68" 11 "82") cols(4))
