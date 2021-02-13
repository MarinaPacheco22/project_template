

mkdir software cod/lib
chmod 755 cod/lib


#Librerias R

PackCran = "'scatterplot3d' 'foreign' 'FSelector' 'MASS' 'GA' 'mclust' 'Hmisc'"

for PackCran in $PackCran
    do 
        Rscript -e 'install.packages('$PackCran', repos="https://cran.rstudio.com/", lib="cod/lib")'
        Rscript -e 'library('$PackCran', lib.loc="cod/lib")'
    done
