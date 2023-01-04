# DiffFracSeq: A Bayesian Model for the Detection of Differential Fractionation of Sequencing Data

DiffFracSeq is a Bayesian statistical model specifically designed to detect differential fractionation.
It represents an alternative way of normalising RNA-Seq data sets from different fractions using the additional information that can be gathered when measuring subsets of a complete sample.
The normalising method uses the transcript counts of samples taken before fractionation to enable reliable inference of RNA-Seq batch specific scale factors within the Bayesian model, rather than relying on a priori estimations.
Below is a figure outlining two experiments that investigate differential fractionation and the ability to detect differential fractionation between DiffFracSeq and DESeq2 is compared.
DiffFracSeq is an open source R software package that will enable even more sensitive comparisons of transcript localisation across conditions and cell types.

![](https://github.com/DimmestP/DiffFracSeq/blob/main/diffFracSeq_diagram.png)

## Installing DiffFracSeq

First install [RStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

```
library(devtools)
devtools::install_github("DimmestP/DiffFracSeq")
```

## Example Usage

```
library(DiffFracSeq)

# load example data
data(simulated_two_fraction_counts)

# check DiffFracSeq model will run
train_DiffFracSeq_model(simulated_two_fraction_counts, iter = 100)
```
