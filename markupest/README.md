Supporting documentation and data for the Stata module markupest - markup estimation with Stata.
The module is available on SSC (ssc install markupest) 

In case of bug reporting, issues, malfunctioning please send an email to:
        gabriele.rovigatti@gmail.com
        
/data:
_ klemscombinedbymeasure.csv - it is the BLS-KLEMS data (combined Sectors and Industry KLEMS Multifactor Productivity Tables) used in Hall (2018) 
_ f_DGP3_replica.dta - Stata 13 dataset with simulated data using ACF (2015) DGP 3
_ data_hall18.dta - Stata 13 dataset with BLS-KLEMS data used to replicate analyses similar to Hall (2018)

/dofile:
_ hall_roeger.do - dofile to run Hall and Roeger estimation on BLS-KLEMS data per NAICS code
_ dlw_simu.do - markup estimation à la De Loecker and Warzynski (2012) on simulated data + graph

