# Retrieve sequences from genbank based on a dataset of accession numbers

Downloads sequences from genbank (nucleotide database) for particular
taxa and genes into a folder called `"0.Sequences"`.

## Usage

``` r
sq.retrieve.indirect(acc.table, download.sqs = FALSE)
```

## Arguments

- acc.table:

  An accession table, ideally generated using `acc.table.retrieve`. The
  data.frame must have the Species, Acc, and gene column names.

- download.sqs:

  Logical indicating whether sequences should be downloaded locally or
  returned as a list.

## Value

None

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
 speciesLevel=TRUE
)

sq.retrieve.indirect(test)

} # }
```
