mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/raw

# from my Desktop:
# scp -r ChipseqARD01-* daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/raw/.

# pulling lanes together
 
for f in `ls raw/Chipseq*/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> ~/temp_NOBACKUP/$label.merged.fastq
done

# number of reads: 

for name in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
   echo -en "$name\n"
   cat ~/temp_NOBACKUP/$name.merged.fastq | wc -l
done

# 122006360/4=30501590
# 162299576/4=40574894
# 211203700/4=52800925
# 192335536/4=48083884

screen -R mapping
# do not have gz files - cat on files 
# will map to hg19 as all Fantom data is for hg19 + motifs
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for name in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
  echo -en $name
  cat ~/temp_NOBACKUP/$name.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > ~/temp_NOBACKUP/$name.sam
done

for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
  echo $sample
  samtools view -Sb ~/temp_NOBACKUP/${sample}.sam > ~/temp_NOBACKUP/${sample}_nonSorted.bam
  samtools sort ~/temp_NOBACKUP/${sample}_nonSorted.bam > ~/temp_NOBACKUP/${sample}.bam
  samtools index ~/temp_NOBACKUP/${sample}.bam
done

# after everything is done
for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
rm ~/temp_NOBACKUP/${sample}.sam ~/temp_NOBACKUP/${sample}_nonSorted.bam
done

for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
mv ~/temp_NOBACKUP/${sample}.bam  ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam/.
mv ~/temp_NOBACKUP/${sample}.bam.bai  ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam/.
done 

# making bw files
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
   echo $sample
  EXTEND=140
  CHROMSIZE=/home/daria/data/data_genomes/hg19/chrom.sizes
  # Number of reads
  librarySize=$(samtools idxstats bam/${sample}.bam | awk '{total+=$3} END{print total}')
  echo $librarySize
  # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
  bamToBed -i bam/${sample}.bam | awk -vCHROM=$CHROMSIZE -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while((getline<CHROM)>0){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g $CHROMSIZE -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > ~/temp_NOBACKUP/${sample}.density.gz	
  # Create WIG file
  gunzip -c ~/temp_NOBACKUP/${sample}.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > ~/temp_NOBACKUP/${sample}.wig.gz
  # Create BigWig file, takes quite some time and memory-consuming
  wigToBigWig ~/temp_NOBACKUP/${sample}.wig.gz $CHROMSIZE bw/${sample}.bw
  # Remove intermediate file
  #rm ${sample}.wig.gz ${sample}.density.gz
done

######################
# putting files in respective folders
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bw


###################### calling peaks
# macs
# following this: https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2
screen -R macs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

for pair in T-KANSL3_T-INPUT K-KANSL3_K-INPUT; do
   exp=$(echo $pair | awk '{split($1,a,"_");print a[1]}')
   input=$(echo $pair | awk '{split($1,a,"_");print a[2]}')
   echo -en "$exp\t$input"
   macs2 callpeak -t bam/$exp.bam -c bam/$input.bam -f BAM -n ${exp}_${input} --outdir macs2 --call-summits -g hs --keep-dup 1
done

# filtering dup and making bed files
for sample in T-INPUT K-INPUT; do
   macs2 filterdup -i bam/${sample}.bam --keep-dup=1 -o macs2/${sample}_filterdup.bed
done&

# random extsize
# making pileup in macs2
# shoudl use 140 extension
#macs2 pileup -i macs2/K-KANSL3_filterdup.bed -o macs2/K-KANSL3_filterdup.pileup.bdg --extsize 130 &
#macs2 pileup -i macs2/T-KANSL3_filterdup.bed -o macs2/T-KANSL3_filterdup.pileup.bdg --extsize 130 &

# bedGraphtoBigWig
# claims that chrM is bigger
# removed chrM
CHROMSIZE=chrom.sizes # problems with chromoc=somal sizes
bedGraphToBigWig macs2/K-KANSL3_filterdup.pileup.bdg $CHROMSIZE bw/K-KANSL3_filterdup.pileup.bw
bedGraphToBigWig macs2/T-KANSL3_filterdup.pileup.bdg $CHROMSIZE bw/T-KANSL3_filterdup.pileup.bw

#

# adding the first line to narrowPeak
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
(echo -en "track type=narrowPeak visibility=3 db=hg19 name=\"K-KANSL3_K-INPUT_narrow_peak\" description=\"K-KANSL3_K-INPUT_narrow_peak\"\n" ; cat macs2/K-KANSL3_K-INPUT_peaks.narrowPeak ) > macs2/K-KANSL3_K-INPUT_peaks.browser.narrowPeak
(echo -en "track type=narrowPeak visibility=3 db=hg19 name=\"T-KANSL3_T-INPUT_narrow_peak\" description=\"T-KANSL3_T-INPUT_narrow_peak\"\n" ; cat macs2/T-KANSL3_T-INPUT_peaks.narrowPeak ) > macs2/T-KANSL3_T-INPUT_peaks.browser.narrowPeak

# loading peaks and bw to the browser
cd ~/web/projects/Alex
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/K-KANSL3_K-INPUT_peaks.browser.narrowPeak K-KANSL3_K-INPUT_peaks.browser.narrowPeak
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/T-KANSL3_T-INPUT_peaks.browser.narrowPeak T-KANSL3_T-INPUT_peaks.browser.narrowPeak
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bw/*.bw .

# links
https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-KANSL3_K-INPUT_peaks.browser.narrowPeak
https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-KANSL3_T-INPUT_peaks.browser.narrowPeak

https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/diff_c1_vs_c2_c3.0_cond2.browser.bed

# tarcks
track type=bigWig name="K-KANSL3.bw" description="K-KANSL3" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-KANSL3.bw
track type=bigWig name="T-KANSL3.bw" description="T-KANSL3" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-KANSL3.bw
track type=bigWig name="T-INPUT.bw" description="T-INPUT" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-INPUT.bw
track type=bigWig name="K-INPUT.bw" description="K-INPUT" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-INPUT.bw

##### differential peaks with macs2: https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/differential

macs2 callpeak -B -t bam/K-KANSL3.bam -c bam/K-INPUT.bam -n K-ChIP --outdir macs2/differential --nomodel --extsize 200 &
macs2 callpeak -B -t bam/T-KANSL3.bam -c bam/T-INPUT.bam -n T-ChIP --outdir macs2/differential --nomodel --extsize 200 &

# differential peak calling
# need to know number of tags after filtering
# total tags in treatment after filtering: 32884992; in input: 39818748
# total tags in control after filtering: in control 43796969

tags1=39818748 # for K
tags2=43796969 # for T
folder=macs2/differential
# -g and -l parameters left unchanged
macs2 bdgdiff --t1 $folder/K-ChIP_treat_pileup.bdg --c1 $folder/K-ChIP_control_lambda.bdg --t2 $folder/T-ChIP_treat_pileup.bdg\
   --c2 $folder/T-ChIP_control_lambda.bdg --d1 $tags1 --d2 $tags2 -g 60 -l 120 --o-prefix $folder/diff_c1_vs_c2

######### calling peaks with macs2, no model
# see above for differential peaks --> 68 peaks  

######################################### peakzilla
mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/peakzilla
peakzilla=~/utils/peakzilla/peakzilla-master/peakzilla.py

# bed files for peakzilla were generated with macs2
#for sample in T-INPUT K-INPUT; do
#   macs2 filterdup -i bam/${sample}.bam --keep-dup=1 -o macs2/${sample}_filterdup.bed
#done&

$peakzilla macs2/T-KANSL3_filterdup.bed macs2/T-INPUT_filterdup.bed > peakzilla/T-KANSL3_peakzilla.tsv
$peakzilla macs2/K-KANSL3_filterdup.bed macs2/K-INPUT_filterdup.bed > peakzilla/K-KANSL3_peakzilla.tsv

$peakzilla macs2/T-KANSL3_filterdup.bed > peakzilla/T-KANSL3_peakzilla.noinput.tsv
$peakzilla macs2/K-KANSL3_filterdup.bed > peakzilla/K-KANSL3_peakzilla.noinput.tsv

###################### genomic distribution
#wget ftp://ftp.ensembl.org/pub/release-71/gtf/homo_sapiens/Homo_sapiens.GRCh37.71.gtf.gz
#gunzip -c Homo_sapiens.GRCh37.71.gtf.gz | awk -vOFS='\t' '($2=="protein_coding"&&(($1>=1&&$1<=22)||$1=="Y"||$1=="X")){print "chr"$0}' | sort -k1,1 -k4,4n > Homo_sapiens.GRCh37.71.gtf
#gunzip -c Homo_sapiens.GRCh37.71.gtf.gz | awk -vFS='\t' -vOFS='\t' -vC="/work2/gschub/anais/data/chrom/hg19.chrom.sizes.good" 'BEGIN{while(getline<C){S[$1]=$2}}{if("chr"$1 in S && $2=="protein_coding"){gid=$9;sub(/gene_id "/,"",gid);sub(/";.*/,"",gid);tid=$9;sub(/.*transcript_id "/,"",tid);sub(/";.*/,"",tid);gn=$9;sub(/.*gene_name "/,"",gn);sub(/";.*/,"",gn);print "chr"$1,$4,$5,$7,$3,gid,tid,gn}}' | sort -k1,1 -k2,2n | gzip > hg19_annotation.txt.gz

# annotations are in ~/data/data_genomes/annotations/raw/
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/data/data_genomes/annotations/raw/hg19_genomic_features.txt

bedtools intersect -wb -a <(awk -vOFS="\t" '{print $1,$2+int(($3-$2)/2),$2+int(($3-$2)/2)}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '{count[$7]++}END{for (i in count)print i,"\t"count[i]}'
# for the whole genome
awk '{count[$4]++}END{for (i in count)print i,"\t"count[i]}' $fileB

# overlaping with CPG islands
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/data/data_genomes/annotations/hg19_cgi.txt

bedtools intersect -c -a <(awk -vOFS="\t" '{print $1,$2+int(($3-$2)/2),$2+int(($3-$2)/2)}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '($4==1){count++}END{print count}'

# overlapping with FANTOM peaks

# preprocessing of FANTOM peaks
#(echo -en "chr\tstart\tend\tstrand\tgene_FANTOM\tglobal_promoter_rank\tCAGE_fresh\tCAGE_revived\tCAGE_thawed\n"
#awk -vOFS="\t" '(NR>1){gsub(/\../,"=",$1) ;a=substr($1,0,length($1)-2); split($1,strand,","); split(a,chr,":"); split(chr[2],pos,"="); #if($4!="NA"){split($2,temp,"@");split(temp[2],g,",");gene=g[1];rank=temp[1]}else{gene="NA";rank="NA"} ; print #chr[1],pos[1],pos[2],strand[2],gene,rank,$7,$8,$9}' hg19.cage_peak_tpm_ann_decoded.osc.txt.gz.extract_THP1.tsv ) > FANTOM_table_THP1.txt

# the same for THP1 but no mean, taking only CAGE_fresh into account
cd /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/CRISPR_CAGE_files/data
awk -vOFS="\t" '(NR>1){print $0}' FANTOM_table_THP1.txt | sort -k5,5 -k7,7gr | awk '{if($5!="NA"){if($7!=0){if(NR==1 || $5!=gene){n=1;gene=$5;print $0, "pRank"n}else{n++; print $0,"pRank"n}}else{print $0,0}}else{print $0,"NA"}}' > ~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt

#bedintersect somehow is not so good
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt
bedtools intersect -loj -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '($4=="."){print $0}' | wc -l 

# closest FANTOM5 peaks --> best solution, better than intersect
# if overlaps several --> will report all
bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '($8=="NA")'

bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '{print $8}' | uniq > /tmp/genes.txt

# overlap of genes with RNA-seq
#files are in ../../RNA-seq/data_May2017/DESeq2/

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($3<=0.01 || $5<=0.01 || $7<=0.01)' | head

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($3<=0.01 || $5<=0.01 || $7<=0.01){print $1}' > /tmp/overlap.genes.txt

awk -vOFS="\t" '($3<=0.01 || $5<=0.01 || $7<=0.01){print $2,$3,$4,$5,$6,$7,$8}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt > /tmp/diffregulated.txt


# the changes were not so dramatic --> will look at the distribution
awk -vOFS="\t" '{print $2,$3,$4,$5,$6,$7,$8}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($7<=0.01){print $6}' | median.R -i - -q -n

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($7<=0.01){print $6}' | median.R -i - -q -n

#R
#library(ggplot2)
#t<-read.table("/tmp/diffregulated.txt")

#par(mfrow = c(1,3))
#for (n in c(2,4,6)){
#genes <- scan("/tmp/overlap.genes.txt",character(), quote = "")
#boxplot(t[,n],ylim=c(-4,1),frame=FALSE)
#for (i in 1:length(genes)){
#y<-c(t[grep(paste0("^",genes[i],"$"),t[,1]),n])
#points<-c(rep(1,t[grep(paste0("^",genes[i],"$"),t[,1]),n]))
#lines(jitter(points),y,type="p")
#text(t[grep(paste0("^",genes[i],"$"),t[,1]),n], labels=genes[i], cex= 0.7, pos=4,xpd=NA)
#}
#}

#points<-c(rep(1,length(log(dist.Pc.min[,2])/log(10))))
#lines(jitter(points),y,type="p")

###################### De novo peak analysis
# de novo against all CpG islands with DREME

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
mkdir -p dreme
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
g=/home/daria/data/data_genomes/hg19/hg19.fa

# see different parameters in help: bedtools getfasta -h 
bedtools getfasta -fi $g -s -bed <(awk '{print $1"\t"$2"\t"$3"\t"NR"\t""+"}' $fileA) -fo dreme/KANSL3_only_peaks.fa -name
# negative regions: some CAGE peaks of expressed genes (CAGEfresh>1)
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt
bedtools intersect -v -a <(awk -vOFS="\t" '($7>=1){print $0}' $fileB | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n)  |  awk '{print $1,$2,$3}' | shuf -n 1000 > dreme/negative.txt
# CAGE peaks are tiny, so will add +/- 100 bp
WIN=100
bedtools getfasta -fi $g -s -bed <(awk -vWIN=$WIN '{print $1"\t"$2-WIN"\t"$3+WIN"\t"NR"\t""+"}' dreme/negative.txt) -fo dreme/negative_win$WIN.fa -name

folder=dreme/KANSL3_v1; mkdir -p $folder
dreme -png -oc $folder -p dreme/KANSL3_only_peaks.fa -n dreme/negative_win100.fa -e 0.05

######################## data from Asifa #################
track type=bigWig name="KANSL3_rep2_mES.bw" description="KANSL3 rep2 mES" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251952_NSL3B_mlES_To1x.bw 
track type=bigWig name="KANSL3_rep1_mES.bw" description="KANSL3 rep1 mES" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251951_NSL3A_mlES_To1x.bw
track type=bigWig name="KANSL3_rep1_NPC.bw" description="KANSL3 rep1 NPC" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251959_NSL3A_mlNPC_To1x.bw
track type=bigWig name="KANSL3_rep2_NPC.bw" description="KANSL3 rep2 NPC" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251960_NSL3B_mlNPC_To1x.bw
 
####### getting fastq data
cd ~/data/published_data/mES/Asifa_KANSL3_mES
mkdir -p fastq

#SRR1258425 KANSL3_mES_rep1
#SRR1258426 KANSL3_mES_rep2
#SRR1258415	input_mES_rep1
#SRR1258416	input_mES_rep2
#SRR1258427 input_NPC_rep1
#SRR1258428 input_NPC_rep2
#SRR1258433 KANSL3_NPC_rep1
#SRR1258434 KANSL3_NPC_rep2

# regular prefetch did nto work
/home/sam/software/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch.2.8.2 -v SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434

screen -R mapping2
fastqdump=/home/sam/software/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump.2.8.2
cd ~/data/published_data/mES/Asifa_KANSL3_mES
for SRR in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434 ; do
   echo -en "$SRR\n"
$fastqdump --outdir fastq --split-files /home/daria/ncbi/public/sra/${SRR}.sra
done

# mapping bowtie PE data
# --strata  do not apply in paired-end mode

cd ~/data/published_data/mES/Asifa_KANSL3_mES
mm9=/k/genomes/mm9/index/bowtie_canonical/mm9
for SRR in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434 ; do
   echo -en "$SRR\n"
   bowtie -p 5 -q -m 1 -v 3 --sam --quiet $mm9 -1 fastq/${SRR}_1.fastq -2 fastq/${SRR}_2.fastq > ~/temp_NOBACKUP/${SRR}.sam
done

# had to move the files to clean Tycho
# new location: ~/data/temp_NOBACKUP/

screen -R mapping  #31561.mapping2
cd ~/data/published_data/mES/Asifa_KANSL3_mES
for sample in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434; do
  echo $sample
  samtools view -Sb ~/data/temp_NOBACKUP/${sample}.sam > ~/data/temp_NOBACKUP/${sample}_nonSorted.bam
  samtools sort ~/data/temp_NOBACKUP/${sample}_nonSorted.bam > ~/data/temp_NOBACKUP/${sample}.bam
  samtools index ~/data/temp_NOBACKUP/${sample}.bam
done

######################## motif analysis with uniprobe hg19

# will have to overlap $peaks with all mapped motifs
#ls ~/data/data_genomes/motifs/uniprobe/mast_p03/hg19/NAME.txt.gz
# will count only once even if 2 overlaps
# no sorting since all files are sorted

screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
peaks=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
folderMotifs=~/data/data_genomes/motifs/uniprobe/mast_p03/hg19

touch motifs/counts.differential.txt
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  echo $name
  number=$(bedtools intersect -u -a $peaks -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/counts.differential.txt
done

# do the same for negatives used for dreme de novo


screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

# take the same number of negatives as positives (seems like A does it this way)
shuf -n 156 dreme/negative.txt > dreme/negative156.txt
negative=dreme/negative156.txt
folderMotifs=~/data/data_genomes/motifs/uniprobe/mast_p03/hg19

rm motifs/counts.negative.txt 
touch motifs/counts.negative.txt

# negatives are expanded to be 500 bp in length
# take the same number of negatives as positives (seems like A does it this way)
# good idea to do the same for positives - expand to a certain length [not done yet]
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  echo $name
  number=$(bedtools intersect -u -a <(awk -vOFS="\t" '($3-$2<=500){print $1,$2-int((500-($3-$2))/2),$3+int((500-($3-$2))/2),"NAME","1"}' $negative | sort -k1,1 -k2,2n) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/counts.negative.txt
done

# calculating enrichments and p-values

totalDiff=$(awk 'END{print NR}' macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed)
totalNeg=$(awk 'END{print NR}' dreme/negative156.txt)

awk -vOFS="\t" -vFile="motifs/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],$2,D,N}' motifs/counts.negative.txt | hyper.R -i - > motifs/counts.hyper.txt


#################################





















