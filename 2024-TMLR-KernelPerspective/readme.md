# Publication Details

E. Gonzalez, M. Abudia, M. T. Jury, R. Kamalapurkar, and J. A. Rosenfeld [**The kernel perspective on dynamic mode decomposition,**](https://arxiv.org/abs/2106.00106) under review at TMLR

Software: MATLAB R2020b or newer (needs the `pagemtimes` function introduced in R2020b)

External Dependencies:
1) The cylinder flow example needs [data](http://dmdbook.com/DATA.zip) from the [DMD book by Kutz et al.](http://www.dmdbook.com/) Unzip `DATA.zip` and copy the `FLUIDS` directory to the placeholder `DATA` directory of this repository. The script `Examples\CylinderFlow\cylinderFlowKoopmanDMD.m` will then reproduce the results in the paper.
2) The turbulent flow example needs the file named [`7999.am`](https://libdrive.ethz.ch/index.php/s/lv7dV40oYlkWJiC/download?path=%2F&files=7999.am), from [this fluid flow data set for machine learning](https://doi.org/10.3929/ethz-b-000515488). Download the file `7999.am` and copy it to the placeholder `DATA` directory. The script `Examples\TurbulentFlow\turbulentFlowKoopmanDMD.m` will then reproduce the results in the paper.

Simulation code is made available as-is without any guarantees whatsoever. Every effort has been made to ensure that the code uploaded here reproduces the results presented in the linked publications; however, the authors do not make any guarantees regarding exact reproduction of the said results. If you notice any mismatch, please contact the owner of the repository.
