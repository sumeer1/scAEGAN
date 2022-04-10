# Package names
# Package names

pkgLoad <- function( packages = "favourites" ) {

    if( length( packages ) == 1L && packages == "favourites" ) {
        packages <- c( "ggplot2", "Seurat", "tidyverse", "cowplot", "pdfCluster", "pheatmap", "cluster", "Rtsne", "umap", "e1071"
        )
    }

    packagecheck <- match( packages, utils::installed.packages()[,1] )

    packagestoinstall <- packages[ is.na( packagecheck ) ]

    if( length( packagestoinstall ) > 0L ) {
        utils::install.packages( packagestoinstall,
                             repos = "http://cran.us.r-project.org"
        )
    } else {
        print( "All requested packages from CRAN already installed" )
    }
    if( length( packagestoinstall) > 0L ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
	BiocManager::install(packagestoinstall_2 )	
    } else {
        print( "All requested packages from BIOCONDUCTOR already installed" )
    }

    for( package in packages ) {
        suppressPackageStartupMessages(
            library( package, character.only = TRUE, quietly = TRUE )
        )
    }
}
pkgLoad("favourites")
