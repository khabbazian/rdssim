### An R package for generating and analyzing respondent-driven sampling referral trees/chains
This package was developed to carry out the simulations section of 
``Khabbazian, Mohammad, et al. Novel sampling design for respondent-driven sampling. arXiv preprint arXiv:1606.00387 (2016).''

#### [rdssim Reference manual](http://homepages.cae.wisc.edu/~khabbazian/pdfs/rdssim.pdf)

### Install using the devtools package.
It requires the version of cpp compiler supports c++11 standards, for instance, GNU GCC (>=4.7).
```
install.packages("devtools")
require(devtools)
install_github("khabbazian/rdssim")
require(rdssim)
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
