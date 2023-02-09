## Reference Database Files for Immuno-Whole Exome Sequence (iWES) Genotyping
A manuscript describing this method of genotyping is in preparation, but you may view the workflow at the GitHub repo [dholab/nf-iWES](https://github.com/dholab/nf-iWES).

Immuno-WES genotyping relies on reference allele sequences from the [Immuno-Polymorphism Database (IPD)](https://www.ebi.ac.uk/ipd/mhc/), which are downloaded as part of the IPD-ref-generator workflow. The workflow then produces the following date-stamped reference databases for Rhesus macaque alleles and Cynomolgus/crab-eating macaque alleles:

- `*.gdna.fasta` - These FASTA files contain only those alleles from IPD that have full-length, genomic DNA sequences, i.e., they contain exons _and_ introns.
- `*.cdna.fasta` - These FASTA files contain all the alleles that lack introns, i.e., do not have full-length genomic DNA sequences published on IPD.
- `*.exon2.fasta` - In these FASTA files, exon-only alleles are trimmed down to exon 2 only. This is important for iWES genotyping because it provides a fall-back for genotyping with reference alleles that are not published with introns.
- `*.immunowes.fasta` - The most recent FASTA file with this extension should be used for iWES Genotyping. It contains all alleles from IPD that are full-length genomic DNAs, and for those that are not, it contains only their exon 2 sequences as a fall back. To make the immunowes FASTA, we simply concatenate the `.exon2.fasta` and `.gdna.fasta`.
