# Retrieve sequences from genbank

Downloads sequences from genbank (nucleotide database) for particular
taxa and genes into a folder called `"0.Sequences"`.

## Usage

``` r
sq.retrieve.direct(
  clades = NULL,
  species = NULL,
  genes = NULL,
  db = "itis",
  maxseqs = 1,
  maxlength = 5000
)
```

## Arguments

- clades:

  A vector listing taxonomic groups of interest (character).

- species:

  A vector listing additional species interest (character). This
  argument can be used to define additional target species in the
  ingroup or species to be sampled in the outgroup (character).

- genes:

  A vector listing gene names of interest (character).

- db:

  Follows `db` in `taxize::downsteam`. Choose from `itis`, `gbif`,
  `ncbi`, `worms`, or `bold`.

- maxseqs:

  Maximum number of sequences to retrieve per search (taxa + gene)
  (numeric).

- maxlength:

  Maximum length of the gene sequence (numeric).

## Value

None

## Examples

``` r
if (FALSE) { # \dontrun{
sq.retrieve.direct(
  clades = c("Felis", "Vulpes", "Phoca"),
  species = "Manis_pentadactyla",
  genes = c("ADORA3", "CYTB")
)
} # }
```
