# Frequently Asked Questions (FAQs)

## What settings are used for the BLAST search?

- `blastp`
- Database: `swissprot`
- Input: `fasta_file.seq`
- Alignments: `500`
- Word size: `6`

## What settings are used for COBALT to create the MSA (Multiple Sequence Alignment)?

- Output format: `mfasta`
- End gap open penalty: `5`
- End gap extend penalty: `1`
- Gap open penalty: `11`
- Gap extend penalty: `1`
- BLAST e-value: `0.003`
- No RPS: `T`
- Tree method: `clust`

## Can I use MoNvIso without an internet connection?

No, MoNvIso requires an internet connection to access necessary databases like UniProt and PDB.

## What species' proteins can be modeled?

Currently, MoNvIso only models proteins belonging to *Homo sapiens* as extracted from the UniProt database.

## How can I cite MoNvIso?

Please cite MoNvIso as follows:

Oliva, F., Musiani, F., Giorgetti, A., De Rubeis, S., Sorokina, O., Armstrong, D. J., Carloni, P., & Ruggerone, P. (2023). Modelling eNvironment for Isoforms (MoNvIso): A general platform to predict structural determinants of protein isoforms in genetic diseases. Frontiers in chemistry, 10, 1059593. [https://doi.org/10.3389/fchem.2022.1059593](https://doi.org/10.3389/fchem.2022.1059593)
