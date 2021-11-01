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
Matteo Saponati, Martin Vinck (2021). **Anticipation of spike patterns at the single neuron level and emergence of spike-timing-dependent plasticity** <br/>
doi: link

```
@article{saponati2021anticipation,
  title={Anticipation of spike patterns at the single neuron level and emergence of spike-timing-dependent plasticity},
  author={Saponati, Matteo and Vinck, Martin},
  journal={arXiv preprint arXiv:lalala},
  year={2021}
}
```


