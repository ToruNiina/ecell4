$ErrorActionPreference = "Stop"
conda build conda.recipe
# Copy the conda build to the home dir, such that it can be registerd as an artifact
#ls C:\Miniconda36-x64\conda-bld
mv C:\Miniconda$env:CONDA\conda-bld .
pwd
ls .
