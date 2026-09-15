# Retrieve accession numbers and titles for a genbank search

Retrieve accession numbers and titles for searches in genbank. This
function is useful for exploring the gene-, population-, and
species-level sampling in genbank.

## Usage

``` r
acc.retrieve(
  organism,
  acc.num = FALSE,
  gene = NULL,
  speciesLevel = FALSE,
  npar = 2,
  nSearchesBatch = 499
)
```

## Arguments

- organism:

  The name of a single taxon (character).

- acc.num:

  Logical indicating whether the `organism` is actually a vector or
  accession numbers.

- gene:

  The name of a single gene region (character; optional).

- speciesLevel:

  Whether the result should be a species-level dataset (logical).

- npar:

  Number of parallel searches (the default is probably the best option).

- nSearchesBatch:

  Number of searches per batch

## Value

This function returns an object of class `data.frame` that includes the
following columns. First, `Species` with the species name listed in
GenBank. Second, `Ti` including the title of the accession in GenBank.
Third, `Acc` listing the accession number. Fourth, `gene` including the
name of the target gene.

## Examples

``` r
if (FALSE) { # \dontrun{
acc.retrieve(organism="Vulpes", gene = 'cytb', species=TRUE)
} # }
```
