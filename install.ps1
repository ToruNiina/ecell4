$ErrorActionPreference = "Stop"
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install -q conda-build anaconda-client
