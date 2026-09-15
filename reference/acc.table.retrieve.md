# Retrieve accession numbers and titles for a given combination of species and genes in genbank

Retrieve accession numbers and titles for complex searches conducted in
genbank. Note that this function depends on `acc.retrieve` as it
actually uses `expand.grid` to find all the relevant organism/gene
combinations in genbank.

## Usage

``` r
acc.table.retrieve(
  clades = NULL,
  species = NULL,
  genes = NULL,
  speciesLevel = NULL,
  npar = 2,
  nSearchesBatch = 499
)
```

## Arguments

- clades:

  A vector of clade names (character). Note that this can either be the
  name of a clade (e.g. Apis) or the code in NCBI (e.g. txid7459).

- species:

  A vector of species names (logical). Note that this can either be the
  name of a clade (e.g. Apis) or the code in NCBI (e.g. txid7459).

- genes:

  A vector of gene names (character; optional).

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
acc.table.retrieve(
 clades  = c('Felis', 'Vulpes', 'Phoca'),
 species = 'Manis_pentadactyla' ,
 genes   = c("A2AB","ADORA3","ADRB2","APOB",
            "APP","ATP7","BCHE","BDNF",
            "BMI1","BRCA1","BRCA2","CNR1",
            "COI","CREM","CYTB","DMP1",
            "EDG1","ENAM","FBN1","GHR",
            "IRBP","ND1","ND2","PLCB4",
            "PNOC","RAG1a","RAG1b","RAG2",
            "TTN","TYR1","VWF"),
 speciesLevel = FALSE
)
} # }
```
