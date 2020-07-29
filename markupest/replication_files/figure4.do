// Hall method //
use "${datadir}/data_hall18", clear

egen nnaics = group(naics)
xtset nnaics year

// time-invariant version with graphical representation of correction
bys naics: markupest mkup_h18_ado, inputs( K L M S E ) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) /*
	*/ go(Y) deltago(deltaY) instruments( deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP) hgraph

gr export "${plotdir}/figure_h18.eps", replace

// time-varying version 
bys naics: markupest mkup_h18_tv, inputs( K L M S E ) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) /*
	*/ go(Y) deltago(deltaY) instruments( deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP) timevarying
	
collapse (mean) mkup_h18_tv _psi_ , by(year)
g tweights = year - 2002
g psi_t = _psi_ * tweights

tsset year
tsline mkup_h18_tv if year >= 1988 & year <= 2017, lw(thick) ylab(1.28(.03)1.38) xlab(1988(5)2018, angle(45)) lc(maroon) xtitle("Year") ytitle("Implied values of {&mu}")

gr export "${plotdir}/figure_h18_panelb.eps", replace
