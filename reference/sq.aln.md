# Align sequences

Perform multiple sequence alignment on all the fasta files saved in a
given folder. Alignment is conducted using
`"BiocManager::install("DECIPHER")"`. Note that this function first
optimizes the direction of the sequences, then aligns using
`"BDECIPHER"`, and finally masks the resulting alignment (optionally but
does it per default). The masking step includes removing common gaps
across all the species and removing highly ambiguous positions. The
resulting aligned sequences are stored to a new folder
`"2.Alignments"`.'

## Usage

``` r
sq.aln(
  folder = "1.CuratedSequences",
  FilePatterns = "renamed",
  sqs.object = NULL,
  mask = TRUE,
  maxFractionGapsSpecies = 0.01,
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
  specific pattern wants to be analyzed.

- sqs.object:

  A list of sequences generated from sq.curate. Only use if you're not
  interested in download sequences locally.

- mask:

  Removes ambiguous sites (Logical, TRUE or FALSE).

- maxFractionGapsSpecies:

  Maximum fraction of gaps per species (when masked)

- ...:

  Arguments passed to `"DECIPHER::AlignSeqs"`.

## Value

This function will return an object of class `list` including the
original and renamed sequence alignments. Optionally, this object will
also include the masked alignments.

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
} # }
```
