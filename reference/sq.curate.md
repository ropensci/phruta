# Curate sequences from genbank

After downloading sequences from genbank, this function curates
sequences based on taxonomic information. Note that this function
provides two summary datasets. First, the accession numbers. Second, the
taxonomic information for each species in the database. The taxonomy
strictly follows the gbif taxonomic backbone. The resulting files are
saved to `"1.CuratedSequences"`. The resulting files also have the most
recent curated taxonomy following the gbif (or selected database)
taxonomic backbone.

## Usage

``` r
sq.curate(
  filterTaxonomicCriteria = NULL,
  mergeGeneFiles = NULL,
  database = "gbif",
  kingdom = NULL,
  folder = "0.Sequences",
  sqs.object = NULL,
  removeOutliers = TRUE,
  minSeqs = 5,
  threshold = 0.05,
  ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
```

## Arguments

- filterTaxonomicCriteria:

  A single string of terms (delimited using "\|") listing all the
  strings that could be used to identify the species that should be in
  the dataset (character).

- mergeGeneFiles:

  A named list, with each element being a character vector indicating
  the names of the files in `"0.Sequences"` that need to be combined
  into a single fasta file. For instance, you can use this argument to
  combine CO1 and COI.

- database:

  A name of a database with taxonomic information. Although 'gbif' is
  faster, it only has information for animals and plants. Other
  databases follow taxize::classification.

- kingdom:

  Optional and only used when database='gbif'. Two possible options:
  "animals" or "plants."

- folder:

  The name of the folder where the original sequences are located
  (character).

- sqs.object:

  A list of sequences generated from `sq.retrieve.indirect`. Only use if
  you're not interested in download sequences locally.

- removeOutliers:

  Whether `odseq:odseq` should be used to remove outliers

- minSeqs:

  minimum number of sequences per locus

- threshold:

  Relative to
  [`odseq::odseq`](https://rdrr.io/pkg/odseq/man/odseq.html). Only
  important if `removeOutliers = TRUE`

- ranks:

  The taxonomic ranks used to examine the taxonomy of the species in the
  `0.Sequences` folder.

## Value

This function will return an object of class `list` with the following
elements. First, the curated sequences with original names. Second, the
curated sequences with species-level names. Third, the accession numbers
table. Fourth, a summary of taxonomic information for all the species
sampled in the files.

## Examples

``` r
if (FALSE) { # \dontrun{
sq.retrieve.direct(
  clades = c("Felis", "Vulpes", "Phoca"),
  species = "Manis_pentadactyla",
  genes = c("ADORA3", "CYTB")
)
sq.curate(
  filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
  database = "gbif", kingdom = "animals",
  folder = "0.Sequences"
)
} # }
```
