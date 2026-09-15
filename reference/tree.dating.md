# Tree dating under treePL or

Performs tree dating under `"treePL"` or `"PATHd-8"` based on secondary
calibrations. Note that `"treePL"` or `"PATHd-8"` must be installed in
your PATH. How to install `"PATHd-8"` in mac
`"https://gist.github.com/cromanpa94/a43bc710a17220f71d796d6590ea7fe4"`
and `"treePL"` can be installed using homebrew (brew install
brewsci/bio/treepl). Thanks to Brian O'Meara and Jonathan Chang,
respectively.

## Usage

``` r
tree.dating(
  taxonomyFolder = "1.CuratedSequences",
  phylogenyFolder = "3.Phylogeny",
  ...
)
```

## Arguments

- taxonomyFolder:

  Name of the folder where `"1.Taxonomy.csv"`, created duing the
  `"sq.curate"` step, is stored (character).

- phylogenyFolder:

  Name of the folder where `"RAxML_bipartitions.phruta"`, created duing
  the `"tree.raxml"` step, is stored (character).

- ...:

  Arguments passed to `"geiger::congruify.phylo"`.

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
tree.dating(
  taxonomyFolder = "1.CuratedSequences",
  phylogenyFolder = "3.Phylogeny", scale = "treePL"
)
} # }
```
