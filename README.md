# Linear E-I circuit

this is a repository for this paper:
<br/> "*Spontaneous variability in gamma dynamics described by a linear harmonic oscillator driven by noise*"<br/>
G Spyropoulos, M Saponati, et al <br/>
(under review, Nature Communications)
[[BiorXiv]](https://www.biorxiv.org/content/10.1101/793729v2.abstract)

## dependencies
the code is written in [Python 3.x]
1. [numpy]
2. [matplotlib]

## usage
`git clone` the repository to your home computer.

## structure  
`python linear_ei_main {} {} {}` - insert the values of beta1 and beta2 and the affine transformation matrix A (as a list of 4 elements `[a1,a2,a3,a4]`)

output: the numerical solution of the AR(2) and of the (E,I) system obtained with the affine transformation 

scripts to replicate the results shown in the related figure
* `supp1_main` 
* `supp2_main`  
* `supp3_main`

## citation and credits
Spyropoulos, G., Saponati, M., Dowdall, J. R., Schölvinck, M. L., Bosman, C. A., Lima, B., ... & Fries, P. (2020). **Spontaneous variability in gamma dynamics described by a linear harmonic oscillator driven by noise**. bioRxiv, 793729. <br/>

```
@article{spyropoulos2020spontaneous,
  title={Spontaneous variability in gamma dynamics described by a linear harmonic oscillator driven by noise},
  author={Spyropoulos, Georgios and Dowdall, Jarrod Robert and Sch{\"o}lvinck, Marieke Louise and Bosman, Conrado Arturo and Lima, Bruss and Peter, Alina and Onorato, Irene and Klon-Lipok, Johanna and Roese, Rasmus and Neuenschwander, Sergio and others},
  journal={bioRxiv},
  pages={793729},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}
```


