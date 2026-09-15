# Shared curation logic for sq.curate

Runs the taxonomy lookup, duplicate removal, outlier detection, and
taxonomic filtering that sq.curate needs regardless of whether the input
sequences come from disk (folder) or from an in-memory sqs.object. Both
entry points in sq.curate build a named list of DNAbin sequences
(fastaSeqs) and hand it here, then take care of their own export step
(writing FASTA/CSV files, or returning a list).

## Usage

``` r
curate.sequences.core(
  fastaSeqs,
  fileNames,
  filterTaxonomicCriteria,
  database,
  kingdom,
  removeOutliers,
  minSeqs,
  threshold,
  ranks
)
```

## Arguments

- fastaSeqs:

  A named list of DNAbin sequences, one element per gene file.

- fileNames:

  The file (or gene) names associated with fastaSeqs, used to track
  which sequences came from which source.

- filterTaxonomicCriteria:

  See `sq.curate`.

- database:

  See `sq.curate`.

- kingdom:

  See `sq.curate`.

- removeOutliers:

  See `sq.curate`.

- minSeqs:

  See `sq.curate`.

- threshold:

  See `sq.curate`.

- ranks:

  See `sq.curate`.

## Value

A list with AccDat, Full_dataset, TableCombined, curatedSeqs, and
toRename.
