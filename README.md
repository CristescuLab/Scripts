# Scripts
Container of useful scripts used in the lab:
- [blast_processing.py](#blastprocessingpy)
- [update_taxonomy.sh](#updatetaxonomysh)

## blast_processing.py
Parse a Blast output on format 6 with the following columns: 'qseqid sseqid pident evalue qcovs qlen length staxid stitle'. It assumes that the stitle contains the lineage or the first two fields are species. If lineage is in stitle, it will assume 7 taxon levels: 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.

```
Usage: blast_processing.py [options] blast_file prefix

If in Compute Canada (Graham/Cedar), use:
module load scipy-stack/2018b

Options:
  -h, --help            show this help message and exit
  -o, --output_filtered
                        Output a TSV with the filtered table [default: False]
  -p PIDENT, --pident=PIDENT
                        Minimum percent identity [default: none]
  -e EVAL, --eval=EVAL  Maximum evalue [default: none]
  -q QCOV, --qcov=QCOV  Minimum query coverage [default: none]
  -Q QLEN, --qlen=QLEN  Minimum query length [default: none]
  -l LENGTH, --length=LENGTH
                        Minimum alignment length [default: none]
  -t TAXLEVEL, --taxlevel=TAXLEVEL
                        Taxonomic level to display [default: species]
  -r MIN_READS, --min_reads=MIN_READS
                        Minimum number of reads to retain group [default: 0]
  -P, --plot            Make a barchart with your group [default: False]
  -a TAX_FOR_PATTERN, --tax_for_pattern=TAX_FOR_PATTERN
                        Parental taxonomic level to subset based on pattern
                        [default: none]
  -b PATTERN, --pattern=PATTERN
                        Pattern to subset the tax_for_pattern with [default:
                        none]
  -s SUFFIX_FOR_PLOT, --suffix_for_plot=SUFFIX_FOR_PLOT
                        Suffix for plot (before extension) [default: none]
  -n NTOP, --ntop=NTOP  Number of hits per query [default: none]
  -c, --use_coi         If no special formating in teh database and using COI
                        [default: False]
```
Let's assumed that you have run the following blast command:

`blastn -db nt -query some_file.fasta -outfmt "6 qseqid sseqid pident evalue qcovs stitle" -out output.blast`

This will generate a TAB-delimited file with the following columns: qseqid, sseqid, pident, evalue, qcovs, qlen, length, staxid,  and stitle. Now, let's imagine that you want to filter your output and generate a species count with only hits that:
1. Have higher than 98% identity
2. Have higher than 98% coverage
3. E-value lower than 1E-10
4. Use only hits where the query lenght is greater than 100 bp
5. That the actual alignment is greater than 60
6. You want to focus your analysis at the genus level

Also, you want the outputs to be prefixed by `test`. You will use this script as follows:
`python blast_processing.py -p 98 -q 98 -e 1E-10 -Q 100 -l 60 -t genus output.blast test`

This will create the following files:
1. test_number_of_reads_in_genus.tsv: This will contain the unique counts of reads per species
2. test_Number_unique_genus_per_read.tsv: How many unique genus (and which ones) you recovered
3. test_List_unique_genus.txt: A list of all unique species in your file that passed the filters.

## update_taxonomy.sh
