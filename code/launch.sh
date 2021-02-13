echo "Ejecucion"
if [-d "cod/lib"];
then
    echo "if"
    Rscript sarsCov.R
else
    echo "else"
    sudo bash setup.sh
    Rscript sarsCov.R
fi
echo "Out"
