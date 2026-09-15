# Run Partitionfinder v.1

This function runs partitionfinder v1 within phruta. For now, all
analyses are based on genes. Please note that you need at least two gene
regions to run partitionfinder.

## Usage

``` r
sq.partitionfinderv1(
  folderAlignments = "2.Alignments",
  FilePatterns = "Masked",
  folderPartitionFinder = "2.1.PartitionFinderv1",
  models = "all",
  run = TRUE
)
```

## Arguments

- folderAlignments:

  Name of the folder where the sequences to align are stored
  (character).

- FilePatterns:

  A string that is common to all the target files in the relevant folder
  (character). Note that this argument can be set to `"NULL"` if no
  specific pattern wants to be analyzed.

- folderPartitionFinder:

  Name of the new folder where the output files are stored (string).

- models:

  Models to run in partitionfinder (string).

- run:

  Run partitionfinder?

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
sq.curate(
  filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
  kingdom = "animals", folder = "0.Sequences"
)
sq.aln(folder = "1.CuratedSequences")
sq.partitionfinderv1(
  folderAlignments = "2.Alignments",
  FilePatterns = "Masked",
  models = "all"
)
} # }
```
