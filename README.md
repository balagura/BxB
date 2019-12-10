# BxB - R package for simulating beam-beam effects in van der Meer scans at LHC

## Simulation
The B*B (B-star-B or BxB) simulation splits the "kicked" bunch into
individual "macro-particles" with weights, and traces them in the
accelerator. Every macro-particle propagates between the interaction points
(IPs) according to the specified betatron phase advances. At every IP it is
influenced by the electromagnetic interaction with the opposite "kicker"
bunch. The original shape of each bunch can be specified as an arbitrary
weighted sum of Gaussians with common mean, round or elliptical in X-Y
projection. The number of IPs is also arbitrary. In the end of the simulation
the modifications of the bunch overlap integrals (ie. of the luminosities) due
to beam-beam are reported.

## Usage
The main function is called "beam_beam". Its arguments can be constructed using helper functions
"kicked", "kickers" and "sim". For details and examples, please, type in R 
```
library(BxB)
?beam_beam
?kicked
?kickers
?sim
```

## Installation
### Direct installation from github:
```
git clone https://github.com/balagura/BxB.git
R CMD build BxB
R CMD INSTALL BxB_<version>.tar.gz
```
or
```
sudo R CMD INSTALL BxB_<version>.tar.gz
```
### Using devtools
If you have not done it yet, install `devtools` pakage
```
(sudo) R
install.packages('devtools')
```
Then, install `BxB`:
```
(sudo) R
devtools::install_github('balagura/BxB')
```
## Other ways to run `BxB`
One can also call `BxB` from python, as a C++ function or as a standalone C++ program with the input 
parameters specified in the configuration file. This is described in https://github.com/balagura/beam-beam-simulation-for-vdM-scans-at-LHC. The underlying `BxB` C++ code is always the same.
