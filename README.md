# Example for Bayesian variable selection book chapter

## Data

The genotyping data are protected. We are therefore using a synthetic dataset emulating the real data. 
The expression and replicated genotyping data are in data/replicated_data.RData.

**IMPORTANT NOTE:** this is a large file which is stored using Git Large File Storage. To clone this file along wi\
th the repository, please install Git LFS, e.g., using Homebrew:

``` bash
brew install git-lfs
```

and then initialise it for your account by using:

``` bash
git lfs install
```
(Alternatively, since this is the only large file, it can be downloaded manually from the Github interface.)

We are then updating the real expression data to simulate genetic associations between the synthetic genotyping data and the transcript levels. 

This last step is done as part of this repos using the R package **echoseq**, which may be installed using the folllowing `devtools` command:

```R
devtools::install_github("hruffieux/echoseq")
```

## Algorithm

The package **atlasqtl** used for the analysis may be installed using
 
```R
devtools::install_github("hruffieux/atlasqtl")
```

## Workflow

The scripts should be executed in the following order:

1. `prepare_data.R` 

2. `atlasqtl_example.R` 

3. `atlasqtl_permutations.R` 

4. `fit_spline.R`

5. `plot_manhattan.R`

## Issues

To report an issue, please use the 
[issue tracker](https://github.com/hruffieux/bayesian_variable_selection_book_chapter/issues) 
at github.com.