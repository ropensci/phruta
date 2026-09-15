# Tree inference under RAxML

Performs tree inference under `"RAxML"` for aligned fasta sequences in a
given folder (default is `"2.Alignments"`). Note that you need at least
two gene regions to run a partitioned analysis.

## Usage

``` r
tree.raxml(
  folder = "2.Alignments",
  FilePatterns = "Masked_",
  raxml_exec = "raxmlHPC",
  Bootstrap = 100,
  outgroup,
  partitioned = FALSE,
  ...
)
```

## Arguments

- folder:

  Name of the folder where the sequences to align are stored
  (character).

- FilePatterns:

  A string that is common to all the target files in the relevant folder
  (character). Note that this argument can be set to `"NULL"` if no
  specific pattern wants to be analized.

- raxml_exec:

  Where to find `"RAxML"` or how to run it from the console? (string).

- Bootstrap:

  Number of bootstrap replicates (numeric).

- outgroup:

  A single string of comma-separated tip labels to be used as outgroup
  in `"RAxML"` See `"RAxML"` documentation for more details (character).

- partitioned:

  Whether analyses should be partitioned by gene (Logical).

- ...:

  Arguments passed to `"ips::raxml"`.

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
tree.raxml(
  folder = "2.Alignments", FilePatterns = "Masked",
  raxml_exec = "raxmlHPC", Bootstrap = 100,
  outgroup = "Manis_pentadactyla"
)
} # }
```
