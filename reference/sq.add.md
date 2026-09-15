# Add local sequences to previously downloaded sequences

This function adds sequences located in a particular folder (default is
`"0.AdditionalSequences"`) to fasta files in another folder (default is
`"0.Sequences"`). Files must be in FASTA format and names of the files
should perfectly match the ones in the previously downloaded folder
(e.g. `"0.Sequences"`). This function creates a third new folder
`"0.0.OriginalDownloaded"` containing the originally downloaded
sequences. The sequences in the `"0.Sequences"` folder get replaced with
the combined ones. Note that new sequences can be added even for just a
fraction of the originally downloaded genes.

## Usage

``` r
sq.add(folderDownloaded = "0.Sequences", folderNew = "0.AdditionalSequences")
```

## Arguments

- folderDownloaded:

  Name of the folder with downloaded sequences.

- folderNew:

  Name of the folder with local sequences.

## Value

None

## Examples

``` r
if (FALSE) { # \dontrun{
sq.add(folderDownloaded = "0.Sequences",
       folderNew = "0.AdditionalSequences")
} # }
```
