## Setting up files for downstream XX-seq analyses

### 1. Retaining only longests CDS per gene
Required files:
> Raw files
>  <ol>
>   <li>faba_reference_qc.R (rename this whatever you want, I started out with faba so here it is)</li>
>   <li>fasta file of the CDS for your organism of interest. Get from Ensembl (i.e. "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/vicia_faba/", then select either CDS or pep .fa file </li>
>  </ol>
Required packages:
>R: 
> <ol>
>   <li>dplyr </li>
>   <li> seqinr </li>
>  </ol>

Download the raw files and run code in RStudio.

