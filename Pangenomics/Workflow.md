### Will house my workflow for analysis of panfaba genome

21/01/26

Ok so first of, we start from (TPC2019-REV-00279R1_Supplemental_Data_Set_1 m (Nod related genes)), based on https://doi.org/10.1105/tpc.19.00279).

Great. This resource offers a great starting point to prominent nodulation-related genes in legumes. However, it's quite the ecletic mixture of genes, with studied genes across various genes being presented in a single column instead of a unified identifier (i.e. Medicago ortholog as identifying column).

What we need to do therefore is to first work out what we need to find matching orthologs for, per gene that has been presented/annotated of interest. 

From https://doi.org/10.1186/1471-2164-13-104, we can work out that Medicago is most closely related to Faba

    "Since M. truncatula is the model legume species that is most closely related to field pea and faba bean"

We therefore need to first get matching Medicago orthologs, which will find the Orthogroups for, and from the Orthogroups then ID what we're interested in for Faba. These genes ID'ed in faba can then be checked across the 51 genomes."

So here's what we do:
1. Get list of Medicago genes from TPC2019-REV-00279R1_Supplemental_Data_Set_1 m (Nod related genes)), based on https://doi.org/10.1105/tpc.19.00279
2. Download the MT4 to MT5 database from https://www.legoo.org/
3. Using this MT4-to-MT5 db, convert names to MT5 so you can look it up in generated OrthoFinder groups.
