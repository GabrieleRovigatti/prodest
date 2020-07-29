/* Graphs on simulated data */
/* In-body + Appendix graphs */
forv DGP = 1(1)3{
	use "${datadir}/f_DGP`DGP'", clear
	tw (kdensity mkupD`DGP'_m1_t if mkupD`DGP'_m1_t < 4, lw(medthick) lc(ebblue)) (kdensity mkupD`DGP'_m2_t if mkupD`DGP'_m2_t < 4, lw(medthick) lp(_) lc(maroon)) /*
		*/ (kdensity mkupD`DGP'_m3_t if mkupD`DGP'_m3_t < 4, lw(medthick) lp(-) lc(forest_green)) (kdensity mkupD`DGP'_m4_t if mkupD`DGP'_m4_t < 4, lw(medthick) lp(dot) lc(sand)) /*
		*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
		*/ xtitle("Markup")
	gr export "${plotdir}/DGP`DGP'_t_d.eps", replace
	
	tw (kdensity mkupD`DGP'_m1_cb if mkupD`DGP'_m1_cb < 4, lw(medthick) lc(ebblue)) (kdensity mkupD`DGP'_m2_cb if mkupD`DGP'_m2_cb < 4, lw(medthick) lp(_) lc(maroon)) /*
		*/ (kdensity mkupD`DGP'_m3_cb if mkupD`DGP'_m3_cb < 4, lw(medthick) lp(-) lc(forest_green)) (kdensity mkupD`DGP'_m4_cb if mkupD`DGP'_m4_cb < 4, lw(medthick) lp(dot) lc(sand)) /*
		*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
		*/ xtitle("Markup") 
	gr export "${plotdir}/DGP`DGP'_cb_d.eps", replace
	
	tw (kdensity mkupD`DGP'_m1_t_cor if mkupD`DGP'_m1_t_cor < 4, lw(medthick) lc(ebblue)) (kdensity mkupD`DGP'_m2_t_cor if mkupD`DGP'_m2_t_cor < 4, lw(medthick) lp(_) lc(maroon)) /*
		*/ (kdensity mkupD`DGP'_m3_t_cor if mkupD`DGP'_m3_t_cor < 4, lw(medthick) lp(-) lc(forest_green)) (kdensity mkupD`DGP'_m4_t_cor if mkupD`DGP'_m4_t_cor < 4, lw(medthick) lp(dot) lc(sand)) /*
		*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
		*/ xtitle("Markup") 
	gr export "${plotdir}/DGP`DGP'_t_cor_d.eps", replace
	
	tw (kdensity mkupD`DGP'_m1_cb_cor if mkupD`DGP'_m1_cb_cor < 4, lw(medthick) lc(ebblue)) (kdensity mkupD`DGP'_m2_cb_cor if mkupD`DGP'_m2_cb_cor < 4, lw(medthick) lp(_) lc(maroon)) /*
		*/ (kdensity mkupD`DGP'_m3_cb_cor if mkupD`DGP'_m3_cb_cor < 4, lw(medthick) lp(-) lc(forest_green)) (kdensity mkupD`DGP'_m4_cb_cor if mkupD`DGP'_m4_cb_cor < 4, lw(medthick) lp(dot) lc(sand)) /*
		*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
		*/ xtitle("Markup") 
	gr export "${plotdir}/DGP`DGP'_cb_cor_d.eps", replace
}
/* helpfile example dataset building and code */
use "${datadir}/f_DGP3", clear
keep if _n < 50001
replace firm = firm + ((rep -1) * 1000)
drop mkup* inv exit lnw lnp rep
xtset firm period
compress
save "${datadir}/f_DGP3_replica", replace
/* run the markup estimation 	*/
forv m = 1/4{
	markupest mkupD3_m`m'_t_corr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm`m') prodestopt("poly(3) acf trans va") corrected verbose
}
/* Plot the graph */
tw (kdensity mkupD3_m1_t_cor if mkupD3_m1_t_cor < 4, lw(medthick) lc(ebblue)) /*
	*/ (kdensity mkupD3_m2_t_cor if mkupD3_m2_t_cor < 4, lw(medthick) lp(_) lc(maroon)) /*
	*/ (kdensity mkupD3_m3_t_cor if mkupD3_m3_t_cor < 4, lw(medthick) lp(-) lc(forest_green)) /*
	*/ (kdensity mkupD3_m4_t_cor if mkupD3_m4_t_cor < 4, lw(medthick) lp(dot) lc(sand)) /*
	*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
	*/ xtitle("Markup") 
