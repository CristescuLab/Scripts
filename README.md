
# Scripts
Container of useful scripts used in the lab:
- [blast_processing.py](#blast_processingpy)
- [submit_blast_n_process.sh](#submit_blast_n_processsh)
- [[launch_submit.sh](#launch_submitsh")
- [update_taxonomy.sh](#update_taxonomysh)

## blast_processing.py
Parse a Blast output on format 6 with the following columns: 'qseqid sseqid pident evalue qcovs qlen length staxid stitle'. It assumes that the stitle contains the lineage or the first two fields are species. If lineage is in stitle, it will assume 7 taxon levels: 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.

```python
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

```bash
blastn -db nt -query some_file.fasta -outfmt "6 qseqid sseqid pident evalue qcovs stitle" -out output.blast
```

This will generate a TAB-delimited file with the following columns: qseqid, sseqid, pident, evalue, qcovs, qlen, length, staxid,  and stitle. Now, let's imagine that you want to filter your output and generate a species count with only hits that:
1. Have higher than 98% identity
2. Have higher than 98% coverage
3. E-value lower than 1E-10
4. Use only hits where the query lenght is greater than 100 bp
5. That the actual alignment is greater than 60
6. You want to focus your analysis at the genus level

Also, you want the outputs to be prefixed by `test`. You will use this script as follows:

```
python blast_processing.py -p 98 -q 98 -e 1E-10 -Q 100 -l 60 -t genus output.blast test
```

This will create the following files:
1. test_number_of_reads_in_genus.tsv: This will contain the unique counts of reads per species
2. test_Number_unique_genus_per_read.tsv: How many unique genus (and which ones) you recovered
3. test_List_unique_genus.txt: A list of all unique species in your file that passed the filters.

## submit_blast_n_process.sh
This is the basic submission for blast and postprocess it. It is supposed to be used with sbatch, providing the filepath where the output files from the QC pipeline are (assume files end in `.trimmed.derep.fasta`), and db_path, the path to the datase (including the database name) you want to use.
Let's assume that your results from the QC pipeline are in this path: `/home/user/work/run`, and that your database of interest is the called `COI` in the path `/home/user/Databases`. Then, you can submit the job my typing:

```bash
sbatch --export=file_path=/home/user/work/run,db_path=/home/user/Databases/COI submit_blast_n_process.sh
```

Giving this, you can submit a large number of jobs that are contained in diferent folders. For example, let's assume that you processed 10 samples, and the outputs are stored in folders called output.some.additional.info.sample1, output.some.additional.info.sample.sample2, etc. You can submit them all with this bash loop:

```bash
for i in `find . -type d -name "output*"`; do
name=`cut -d'.' -f5 ${i}`
sbatch --job-name=${name} --export=file_path=${i},db_path=/home/user/Databases/COI submit_blast_n_process.sh
done
```

## update_taxonomy.shAll this will be done with 8 threads on the blast. If you want to modify this, you can pass it to the export option of sbatch. Imagine you want to use 32 threads, then you'll have to submitted like this:

```bash
sbatch --cpus-per-task=32 --export=file_path=/home/user/work/run,db_path=/home/user/Databases/COI,cpus=32 submit_blast_n_process.sh
```

## launch_submit.sh
This script is a simple for loop to launch multiple scripts  jobs (as in the previous example). It requires that the files to be processed be in a subfolder within the one where you will be launching the job. It assumes that your folder is named with the last dot separated field being your sample name (or the prefix you want your job to run in).

### Usage:  
bash launch_submit.sh <pattern_folder> <database_with_path> <path_to_code> \<cpus> \<mem> \<account>  

\<cpus> \<mem> \<account> are optional, defaulting to 8, 32G and def-mcristes (for Cristescu lab members)

#### Examples:
Lets imagine that you have to launch multiple jobs on the same submission script called `submit.sh`, that is in the path `/home/user/`. Lets also imagine that you have multiple folders (let's say 10) named `output.sample1`, `output.sample2`, ... ,  `output.sample10`, and within them a file (or files) that will be called by  `submit.sh`. Also, say your submitfile has a variable `db_path` that is `/home/user/Databases/NCBI/nt`. For now, let us assume that you want the default threads (8), memory (32G), and account (def-mcristes), then you can submit all folders in one go by:

```bash
bash launch_submit.sh output /home/user/Databases/NCBI/nt /home/user/submit.sh
```

This will submit everyone of your folders to a different job in the cluster. Now, lets say you want to modify the cpus to work with, the memory, and the account. Let's say that you want 32 cpus, on 146G of memory, and using the project rrg-professor. Then you will do:

```bash
bash launch_submit.sh output /home/user/Databases/NCBI/nt /home/user/submit.sh 32 146G rrg-professor
```

## update_taxonomy.sh
This script will update the accession2taxid files, and will create the `accession2taxid.lineage` file with he following structure:

|  |  |
|--|--|
|  |  |

So it maps the accession numbers with the respective taxonomic lineage. **BE AWARE** this scripts takes lots of time, cpu and memory. Use it very sparingly to make sure that the taxonomy remains up to date.

The usage is simple:

`bash update_taxonomy.sh`

You should be in the folder that you want the files to be (and where mostlike they are already). For our lab, this files are in the `Databases` folder, under the subfolder `accession2taxid`

## sync_folder.sh
Simple script to sync two servers on a particular folder. This script only works if you have an account in **BOTH** servers. The sync will update the files that need to be updated leaving the rest untouch.
the usage is:

` bash sync_software.sh <path2folder> <source_server> <target_server>`

So, for example, let's imagine that you want to sync the `Databases` folder that is the the path `~/projects/def-mcristes/Databases/`. Let's assume that you have the `sync_folder.sh` in your home directory `~`, and that you are in the `graha,` cluster and want to sync to the `cedar` cluster. First go to the folder you wan to sync and then execute the file:
```
cd ~
 bash ~/sync_software.sh `~/projects/def-mcristes/Databases/  graham.computecanada.ca cedar.computecanada.ca
```
** You need to type the password for the graham first and then for cedar **
