# RogueNaRok within phruta

Implements the RogueNaRok algorithm for rogue taxon identification
within phruta

## Usage

``` r
tree.roguetaxa(folder = "3.Phylogeny", ...)
```

## Arguments

- folder:

  Name of the folder where the sequences to align are stored
  (character).

- ...:

  Arguments passed to `"Rogue::RogueTaxa"`.

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
tree.roguetaxa(folder = "3.Phylogeny")
} # }
```
