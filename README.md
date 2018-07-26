# proteogenomics_python
python scripts for proteogenomics analysis.

The whole workflow has been automated into one [nextflow](https://github.com/lehtiolab/proteogenomics-analysis-workflow) pipeline.

map_peptide2genome.py is a python script to map known peptides back to genome.You need three input files:  
a gtf annotation file  
a fasta file including protein sequences  
a IDmap file which contains gene id, transcript id and protein id  

IDmap file can be downloaded using Ensembl Biomart tool. See IDmap_file_example.txt. 

`python map_peptide2genome.py --input input_filename --gtf Homo_sapiens.GRCh37.ensembl87.gtf --fasta Homo_sapiens.GRCh37.ensembl87.pep.all.fa  --IDmap Ensembl87_IDlist.txt --output output_filename`

--input: peptide sequence in first column, protein accession in second column

3frame_translation.py is  python script to do three frame translation.(default standard code)
Example:

`python 3frame_translation.py genome.fasta genome.3FT.fasta`

sixframetranslation.py is python script to do six frame translation and full trypsin digestion at the same time.
Example: 

`python sixframetranslation.py --input genome.fasta --output genome.6FT.txt --nuclear_trans_table 1 --mito_trans_table 2 --min_length 8 --max_length 30`

--nuclear_trans_table is to specify translation table used for nuclear DNA
--mito_trans_table is to specify translation table used for mitochondrial DNA



# Curation of novel peptides from VarDB search
1. Map novel peptides back to genome

`python map_novelpeptide2genome.py --input novpep.txt --gtf VarDB.gtf --fasta VarDB.fasta --gff_output example_novpeps.gff3 --tab_out example_novpep.hg19cor.txt`

The input file `novpep.txt` must contain the two columns with the name: `Peptide` and `Protein`, which are required to map them to genome.

2. Make fasta file for novel peptides

`python to_fasta.py example_novpeps.txt example_novpeps.fasta`

3. BLASTP analysis

`blastp -db UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta -query ../example_novpep.fasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 8 -max_target_seqs 1 -evalue 1000 -out example_novpeps.blastp.out.txt`

4. Parse BLASTP output

`python parse_blastp_out.py --input example_novpep.hg19cor.txt --blastp_result example_novpep.blastp.out.txt --fasta UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta --output example_novpeps.blastp.parsed.txt`

5. Annotate loci - annovar

`python prepare_annovar_input.py --input example_novpep.hg19cor.txt --output example_novpep_avinput.txt`

you need to install annovar before you can run the next command. [Annovar Download page](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

`./annotate_variation.pl -out example_novpep -build hg19 example_novpep_avinput.txt humandb/`

6. Parse annovar result

`python parse_annovar_out.py --annovar_out example_novpep.variant_function --input example_novpeps.blastp.parsed.txt --output example_novpeps.blastp.annovar.txt`

7. Extract PSMs of novel peptide with single substitution

The column name of peptide sequence should be "Peptide", otherwise use --peptide_column to specify a different name

`python extract_1mismatch_novelpsm.py example_novpeps.blastp.annovar.txt example_novpeps.psms.txt example_novpep_1mismatch.psm.txt`

8. Run SpectrumAI [Download SpectrumAI here](https://github.com/yafeng/SpectrumAI).

9. Parse SpectrumAI result

`python parse_spectrumAI_out.py --spectrumAI_out specAI_file --input example_novpeps.blastp.annovar.txt --output output_filename`


# Validation in orthognal evidence
1. calculate conservation scores 
    
`python calculate_phastcons.py novel_peptides.gff3 hg19.100way.phastCons.bw novpeps.phastcons.scores.txt`

2. predict phyloCSF coding potential [see description here](https://github.com/hussius/gff-phylocsf-human).

`python calculate_phyloscf.py novel_peptides.gff3 file_path_to_bigwig novpeps.phyloCSF.scores.txt`

3. count reads support for novel peptides in Bam files

`python scam_bams.py --input_gff novel_peptides.gff3 --bam_files bam_files_list.txt --output novelpep_readcount.txt `


## build mutant protein and peptide sequences DB from COSMIC

```
sftp 'your_email_address@example.com'@sftp-cancer.sanger.ac.uk

sftp> get cosmic/grch38/cosmic/v85/CosmicMutantExport.tsv.gz
sftp> get cosmic/grch38/cosmic/v85/All_COSMIC_Genes.fasta.gz
sftp> exit

python convertCOSMIC2_mutant_protein.py All_COSMIC_Genes.fasta CosmicMutantExport.tsv Cosmic_v85_mutant_protein.fasta
python digest_mutant_protein.py All_COSMIC_Genes.fasta Cosmic_v85_mutant_protein.fasta Cosmic_v85_mutant_peptides.fasta
```
