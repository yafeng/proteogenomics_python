# proteogenomics_python
python scripts for proteogenomics analysis

map_peptide2genome.py is a python script to map known peptides back to genome.

python map_peptide2genome.py --input input_filename --gtf Homo_sapiens.GRCh37.75.gtf --fasta Homo_sapiens.GRCh37.75.pep.all.fa  --IDmap Ensembl75_IDlist.txt --output output_filename"

--input: peptide sequence in first column, protein accession in second column

3frame_translation.py is  python script to do 3 frame translation.


#The following scripts are used in the curation stage of the published proteogenomics workflow. 

1. combine vardb and 6FT search result,remove redundant ones

python combine_6RF_var_overlap.py ./example_vardb/peptide_table.txt ./example_6rf/peptide_table.txt example_vardb_6rf_novpep.txt

2. extract PSM table for combined novel peptide table

python extract_6RF_var_overlap_psm.py example_vardb_6rf_novpep.txt ./example_vardb/psm_table.txt ./example_6rf/psm_table.txt example_vardb_6rf_novpep.psm.txt

3. Map novel peptides back to genome

python map_novelpeptide2genome.py --input ../example_vardb_6rf_novpep.txt --gtf VarDB.gtf --fasta VarDB.3frame.fasta --gff_output example_vardb_6rf_novpep.gff3 --tab_output example_vardb_6rf_novpep.hg19cor.txt


4. Make fasta file for novel peptides

python to_fasta.py example_vardb_6rf_novpep.txt example_vardb_6rf_novpep.fasta

5. BLASTP analysis

blastp -db UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta -query ../example_vardb_6rf_novpep.fasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 8 -max_target_seqs 1 -evalue 1000 -out example_vardb_6rf_novpep.blastp.out.txt

6. Parse BLASTP output

python parse_blastp_out.py --input map2genome/example_vardb_6rf_novpep.hg19cor.txt --blastp_result blastpDB/example_vardb_6rf_novpep.blastp.out.txt --fasta blastpDB/UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta --output example_vardb_6rf_novpep.hg19cor.blastp.txt

7. Annotate loci - annovar

python prepare_annovar_input.py --input example_vardb_6rf_novpep.hg19cor.txt --output example_novpep_avinput.txt

you need to install annovar before you can run this following command

./annotate_variation.pl -out example_novpep -build hg19 example_novpep_avinput.txt humandb/

8. Parse annovar result

python parse_annovar_out.py --annovar_out example_novpep.variant_function --input example_vardb_6rf_novpep.hg19cor.blastp.txt --output example_vardb_6rf_novpep.hg19cor.blastp.annovar.txt

9. Extract PSMs of novel peptide with single substitution

The column name of peptide sequence should be "Peptide", otherwise use --peptide_column to specify a different name

python extract_single-sub-novelpep_psm.py --input_psm example_vardb_6rf_novpep.psm.txt --input_pep example_vardb_6rf_novpep.hg19cor.blastp.annovar.txt --peptide_column "Peptide" --output example_novpep_1mismatch.psm.txt

10. Run SpectrumAI


11. Parse SpectrumAI result

python parse_spectrumAI_out.py --spectrumAI_out example_novpep_1mismatch.SpectrumAI.txt --input example_vardb_6rf_novpep.hg19cor.blastp.annovar.txt --output example_vardb_6rf_novpep.hg19cor.blastp.annovar.SpectrumAI.txt

