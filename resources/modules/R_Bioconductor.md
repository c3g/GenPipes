
## Installing and using R ##

### Installation


### Use

To use R, do something like this:

	module load mugqic/R_Bioconductor/3.1.2_3.0
	module load mugqic/mugqic_R_packages/1.0.0
	R --no-restore --no--save
	
R should never be invoked with *--vanilla*. Doing so makes R ignore *etc/Rprofile* which contains a few important setting (e.g bitmapType option)

### For developers


#### Installing a package to mugqic\_R\_packages ####

The paradigm is to install our in-house packages in a different location than CRAN/Bioconductor packages. This simpifies maintaining different versions as development progresses. The modules mugqic\_R\_packages serve the purpose of pre-pending the locaiton of these external libraries to $R\_LIBS.

While developing 





#### Daily deploy of mugqic\_dev/R\_Bioconductor ####



