# SARS-COV-2 RNA editing
<b>A-to-I RNA editing in SARS-COV-2</b>

<p>RNA editing by adenosine deamination is a key molecular mechanism to modulate the innate immune response. During the course of a viral infection, it can have antiviral effects by destabilizing double stranded RNAs (dsRNAs) or proviral effects by suppressing the immune response.<br>
Here we provide the computational workflow to detect A-to-I RNA editing in the SARS-COV-2 genome. Our scripts enable the identification of hyper-edited reads as well as A-to-I editing at single nucleotide resolution.</p><br>
<b>Requirements</b>
<ul>
<li><b>Software</b></li>
<ul>
<li>FASTP for read trimming</li>
<li>BWA for read mapping</li>
<li>GSNAP for read mapping</li>
<li>sambamba for filtering of SAM/BAM files</li>
<li>SAMtools for manipulating BAM files</li>
<li>REDItools v2 for RNA editing calling</li>
<li>bam2fastq to estract reads from BAM files</li>
<li>python 2.x to run scripts and the following mandatory modules:</li>
  <ul>
    <li>Pysam</li>
  </ul>
</ul>
<li><b>Genomic Data</b></li>
<ul>
<li>SARS-COV-2 genome (NC045512.2)</li>
<li>Human Genome hg19</li>
<li>Reads from PRJNA625518 BioProject - RNAseq from infected Calu-3 cells (only paired end runs):</li>
<ul>
<li>4 hours post infection: SRR11550047 - SRR11550048</li>
<li>12 hours post infection: SRR11550043 - SRR11550044</li>
<li>24 hours post infection: SRR11550045 - SRR11550046</li>
</ul>
</ul>
</ul>
<br>
<b>Environment</b>
<p>Create a folder named PRJNA625518 and a subfolder for each run. Put Fastq files in each subfolder</p>
<div class="highlight-python">
<pre>
mkdir PRJNA625518
cd PRJNA625518
mkdir SRR11550046 # repeat for all runs
</pre>
</div>
<p>The folder structure shoudl be:</p>
  <ul>
    <li>PRJNA625518</li>
    <ul>
    <li>SRR11550046</li>
    <ul>
    <li>SRR11550046_1.fastq.gz</li>
    <li>SRR11550046_2.fastq.gz</li>
    </ul>
    </ul>
</ul>
<b>Procedure</b>
<p>Here we illustrate the workflow step by step.</p>
<div class="highlight-python">
<pre>
cd PRJNA625518/SRR11550046
# fastq trimming
fastp -i SRR11550046_1.fastq.gz -o SRR11550046_1_trim.fastq.gz -I SRR11550046_2.fastq.gz -O SRR11550046_2_trim.fastq.gz -w 10 -q 25 -u 20 -l 50 -x --cut_tail --cut_tail_mean_quality 25 --trim_front1 0 --trim_tail1 0
# mapping of reads against hg19+virus
bwa mem -t 8 hg19covid SRR11550046_1_trim.fastq.gz SRR11550046_2_trim.fastq.gz | samtools sort -@ 8 -o SRR11550046_all.bam -
samtools index -@ 8 SRR11550046_all.bam
# extract unique viral reads
sambamba-0.7.1 view -F "not (unmapped or mate_is_unmapped or supplementary or secondary_alignment or chimeric) and proper_pair" SRR11550046_all.bam NC_045512.2 > SRR11550046.sam
cat header SRR11550046.sam | samtools sort -o SRR11550046.bam -
samtools index SRR11550046.bam
rm SRR11550046.sam
# calculate viral coverage depth
samtools depth -a SRR11550046.bam > SRR11550046.dp
python getCoverage.py SRR11550046.dp > SRR11550046.dp.out
# infer strand orientation
infer_experiment.py -i SRR11550046.bam -r covid_plus.bed > SRR11550046.inf
# extract viral reads in fastq format
bam2fastq -o SRR11550046#.fastq SRR11550046.bam
# align viral reads by GSNAP
gsnap -D gmapdb -d sarscov2 -c sarscov2tr -s splicesites.iit --nofails -n 1 -Q -N 1 -A sam -t 8 SRR11550046_1.fastq SRR11550046_2.fastq | samtools sort -@ 8 -o SRR11550046.gsnap_temp.bam
samtools index -@ 8 SRR11550046.gsnap_temp.bam
sambamba-0.7.1 view -t 8 -h -F "not (unmapped or mate_is_unmapped or supplementary or secondary_alignment or chimeric) and proper_pair" SRR11550046.gsnap_temp.bam | samtools sort -@ 8 -o SRR11550046.gsnap.bam
samtools index -@ 8 SRR11550046.gsnap.bam
rm SRR11550046_1.fastq SRR11550046_2.fastq
rm SRR11550046.gsnap_temp.bam*
mv SRR11550046.bam SRR11550046.bwa.bam
mv SRR11550046.bam.bai SRR11550046.bwa.bam.bai
</pre>
</div>
<p>Detect hyper-edited reads</p>
<div class="highlight-python">
<pre>
# use the script SubstitutionsPerSequence.py. It takes 3 parameters:
# BAM file of aligned viral reads
# Info about the strand orientation. 0: unstranded 1: FR Second Strand 2: FR First Strand
# estimate the fraction of read to filter AG cluster y: estimate (for long reads) n: fix it at 0.05
python SubstitutionsPerSequence.py SRR11550046.bwa.bam 2 n
</pre>
</div>
<p>SubstitutionsPerSequence.py generates 6 output files:</p>
<ul>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
</ul>
/opt/exp_soft/biomed/epicardi/python/bin/python /lustrehome/epicardi/home/covid/runREDItools.py Vero_SCV2-0.bam


