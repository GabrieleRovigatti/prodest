use "${datadir}/f_DGP3_replica.dta", clear
/* run the markup estimation 	*/
forv m = 1/4{
	markupest mkupD3_m`m'_t_corr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm`m') prodestopt("poly(3) acf trans va") corrected verbose
}
/* Run the graph */
tw (kdensity mkupD3_m1_t_cor if mkupD3_m1_t_cor < 4, lw(medthick) lc(ebblue)) /*
	*/ (kdensity mkupD3_m2_t_cor if mkupD3_m2_t_cor < 4, lw(medthick) lp(_) lc(maroon)) /*
	*/ (kdensity mkupD3_m3_t_cor if mkupD3_m3_t_cor < 4, lw(medthick) lp(-) lc(forest_green)) /*
	*/ (kdensity mkupD3_m4_t_cor if mkupD3_m4_t_cor < 4, lw(medthick) lp(dot) lc(sand)) /*
	*/ , legend(order( 1 "No Meas Error" 2 "{&sigma}{superscript:2}{subscript:m}=0.1" 3 "{&sigma}{superscript:2}{subscript:m}=0.2" 4 "{&sigma}{superscript:2}{subscript:m}=0.5")) /*
	*/ xtitle("Markup") 
