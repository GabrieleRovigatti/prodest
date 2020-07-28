use "${datadir}/data_hall18", clear

egen nnaics = group(naics)
xtset nnaics year

g naics_2d = substr(naics, 1, 2)	

// time-invariant Hall version 
bys naics_2d: markupest mkup_h18, inputs( K L M S E ) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) /*
	*/ go(Y) deltago(deltaY) instruments( deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP) 
// Roeger
bys naics_2d: markupest mkup_roeg, inputs(L M S E  K) method(roeger) go(Y)
// last 10 years
bys naics_2d: markupest mkup_h18_10y if year > 2007, inputs( K L M S E ) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) /*
	*/ go(Y) deltago(deltaY) instruments( deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP) 
bys naics_2d: markupest mkup_roeg_10y if year > 2007, inputs(L M S E K) method(roeger) go(Y) 
	
collapse (mean) mkup_h* mkup_r*, by(naics_2d)

replace naics_2d = "Agriculture, Forestry, Fishing, Hunting" if naics_2d == "11"
replace naics_2d = "Mining" if naics_2d == "21"
replace naics_2d = "Utilities" if naics_2d == "22"
replace naics_2d = "Construction" if naics_2d == "23"
replace naics_2d = "Manufacturing" if naics_2d == "31"
replace naics_2d = "Manufacturing" if naics_2d == "32"
replace naics_2d = "Manufacturing" if naics_2d == "33"
replace naics_2d = "Wholesale Trade" if naics_2d == "42"
replace naics_2d = "Retail Trade" if naics_2d == "44"
replace naics_2d = "Transportation and Warehousing" if naics_2d == "48"
replace naics_2d = "Transportation and Warehousing" if naics_2d == "49"
replace naics_2d = "Information" if naics_2d == "51"
replace naics_2d = "Finance and Insurance" if naics_2d == "52"
replace naics_2d = "Real Estate Rental and Leasing" if naics_2d == "53"
replace naics_2d = "Professional, Scientific, and Technical Services" if naics_2d == "54"
replace naics_2d = "Management of Companies and Enterprises" if naics_2d == "55"
replace naics_2d = "Administrative and Support and Waste Manag" if naics_2d == "56"
replace naics_2d = "Educational Services" if naics_2d == "61"
replace naics_2d = "Health Care and Social Assistance" if naics_2d == "62"
replace naics_2d = "Arts, Entertainment, and Recreation" if naics_2d == "71"
replace naics_2d = "Accommodation and Food Services" if naics_2d == "72"
replace naics_2d = "Other Services (No Public Admin)" if naics_2d == "81"
 
eststo fhall: estpost tabstat mkup_h18, by(naics_2d) stat(mean) nototal listwise elabels
eststo froeger: estpost tabstat mkup_roeg, by(naics_2d) stat(mean) nototal listwise elabels

eststo fhall_10: estpost tabstat mkup_h18_10, by(naics_2d) stat(mean) nototal listwise elabels
eststo froeger_10: estpost tabstat mkup_roeg_10, by(naics_2d) stat(mean) nototal listwise elabels

esttab fhall froeger fhall_10 froeger_10 using "${tabledir}/table1.tex", f tex main(mean) wide nostar label replace coeflabels(`e(labels)') noobs nonum /*
	*/ mti("$\hat{\mu}^{hall}$" "$\hat{\mu}^{roeger}$" "$\hat{\mu}^{hall}$" "$\hat{\mu}^{roeger}$") /*
	*/ mgroups("1987-2017" "2008-2017" , pattern(1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cline{@span})) 
