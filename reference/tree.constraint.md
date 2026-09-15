# Tree inference under RAxML

Performs tree inference under `"RAxML"` for aligned fasta sequences in a
given folder (default is `"2.Alignments"`).

## Usage

``` r
tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus",
    "species_names"),
  Topology = "((ingroup), outgroup);",
  outgroup = NULL
)
```

## Arguments

- taxonomy_folder:

  Name of the folder where the 1.Taxonomy file is stored.

- targetColumns:

  Where to find `"RAxML"` or how to run it from the console? (string).

- Topology:

  A string summarizing the desired topological constraint in newick
  format.

- outgroup:

  Optional and only required when the topology argument is not
  "((ingroup), outgroup);".

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

tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family",
  "genus", "species_names"),
  Topology = "((ingroup), outgroup);",
  outgroup = "Manis_pentadactyla"
)
tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family",
  "genus", "species_names"),
  Topology = "((Felis), (Phoca));"
)
} # }
```
