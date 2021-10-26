####Trimmomatic cleaning example script ##################

for j in `ls *_1.fq.gz | sed 's/_1.fq.gz//' `
do
  echo $j; read
  nombre1=$j"_1.fq.gz"
  nombre2=$j"_2.fq.gz"
  java -jar /usr/local/bin/trimmomatic-0.36.jar \
  PE -threads 6 ${j}_1.fq.gz ${j}_2.fq.gz \
  scleaned_pe_paired/$nombre1 scleaned_pe_unpaired/$nombre1 \
  scleaned_pe_paired/$nombre2 scleaned_pe_unpaired/$nombre2 \
  HEADCROP:19
done

####Prinseq ducplicate removal example script ############

for j in `ls *_1.fq | sed 's/_1.fq//' `
do
  echo $j; read
  nombre1=$j"_1.fq.gz"
  nombre2=$j"_2.fq.gz"
  perl /tmp/prinseq-lite-0.20.4/prinseq-lite.pl \
  -derep 12345 -fastq ${j}_1.fq -fastq2 ${j}_2.fq
done

###

bowtie2-build CASILVA_DNA.fasta CAT_SILVA_DNA
seqmagick convert a.fasta a.sto


#####samples_file.txt######################################
cond_A    cond_A_rep1	Leaft-1_1.fq	Leaft-1_2.fq    
cond_B    cond_B_rep1	SC1-1_1.fq	SC1-1_2.fq
cond_B    cond_B_rep2	SC2-1_1.fq	SC2-1_2.fq
cond_B    cond_B_rep3	SC4-1_1.fq	SC4-1_2.fq
cond_C    cond_C_rep1	Root-1_1.fq	Root-1_2.fq
cond_D    cond_D_rep1	ST1-1_1.fq	ST1-1_2.fq
cond_D    cond_D_rep2	ST3-1_1.fq	ST3-1_2.fq
cond_D    cond_D_rep3	ST4-1_1.fq	ST4-1_2.fq

#########Triniy de novo assembly example script############
Trinity --seqType fq --max_memory 90G  \
--samples_file samples_file.txt --CPU 7 \
--output /trinity_out


##Alignment of clean reads to assembly example script#######
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta \
-1 ../Leaft-1_1.fq -2 ../Leaft-1_2.fq \ 
\2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam


##BLASTX to Uniprot for qual check #########################
blastx -query /Trinity.fasta -db uniprot_sprot.fasta \
-out blastx.outfmt6 -evalue 1e-20 -num_threads 6 \
-max_target_seqs 1 -outfmt 6

analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 \
Trinity.fasta uniprot_sprot.fasta

blast_outfmt6_group_segments.pl blastx.outfmt6 Trinity.fasta \
uniprot_sprot.fasta >blastx.outfmt6.grouped

blast_outfmt6_group_segments.tophit_coverage.pl \
blastx.outfmt6.grouped


##### BUSCO analysis example scritp #########################
busco -m transcriptome -i Trinity.fasta  -o OUTPUT_busco \
-l eudicots_odb10

busco -m transcriptome -c 12 -i Trinity.fasta  -o OUTPUT_busco \
--lineage_dataset ./eudicots_odb10

######CD-HIT#################################################
cd-hit-est -i Trinity.fasta -o Unigenes_cdhit.fasta -c 0.95 -T 10

######Peptide prediction TRANSDECODER########################
TransDecoder.LongOrfs -t Trinity.fasta --gene_trans_map \
Trinity.fasta.gene_trans_map

######Annotation#############################################
makeblastdb -in plant.207.protein.faa -dbtype prot


sed "s/*//g" Trinity.fasta.transdecoder.pep  \
>Trinity.fasta.transdecoder.asterisc.pep 


blastp -query ../Trinity.fasta.transdecoder.asterisc.pep\
 -db plant.207.protein.faa -out nr_specific_blastx.outfmt6\
 -evalue 1e-5 -num_threads 15 -max_target_seqs 1 \
 -outfmt '6 qaccver saccver pident qcovs length qstart \
 qend sstart send evalue ssciname  stitle'


########Trnascript Abundance estimation RSEM################
align_and_estimate_abundance.pl  --transcripts Trinity.fasta \
--seqType fq --samples_file samples_file.txt  --est_method RSEM\
 --output_dir RSEM  --aln_method bowtie2 --thread_count 3 \
 --gene_trans_map Trinity.fasta.gene_trans_map --prep_reference
