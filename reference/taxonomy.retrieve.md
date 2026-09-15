# Curate sequences from genbank

After downloading sequences from genbank, this function curates
sequences based on taxonomic information. Note that this function
provides two summary datasets. First, the accession numbers. Second, the
taxonomic information for each species in the database. The taxonomy
strictly follows the gbif taxonomic backbone. The resulting files are
saved to `"1.CuratedSequences"`. The resulting files also have the most
recent curated taxonomy following the gbif taxonomic backbone.

## Usage

``` r
taxonomy.retrieve(
  species_names = NULL,
  database = "gbif",
  kingdom = NULL,
  ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
```

## Arguments

- species_names:

  A vector of species names

- database:

  A name of a database with taxonomic information. 'gbif' only works for
  animals and plants. Databases follows taxize::classification

- kingdom:

  Optional and only used when database = 'gbif'. Two possible options:
  "animals" or "plants."

- ranks:

  Taxonomic ranks to retrieve from the target species.

## Value

`data.frame` of taxonomic information for the target species (valid in
the database)

## Examples

``` r
if (FALSE) { # \dontrun{
taxonomy.retrieve(
  species_names = c(
    "Felis_catus", "PREDICTED:_Vulpes",
    "Phoca_largha", "PREDICTED:_Phoca",
    "PREDICTED:_Manis", "Felis_silvestris", "Felis_nigripes"
  ),
  database = "gbif", kingdom = "animals"
)
} # }
```
