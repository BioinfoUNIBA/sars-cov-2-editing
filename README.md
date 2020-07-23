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
<li>SRR11550046.bwa.Badreads with a list of low quality reads</li>
<li>SRR11550046.bwa.context with a list of nucleotides surrounding the editing sites to calculate the sequence context</li>
<div class="highlight-python">
<pre>
#Example of the content
# column 1: genomic position
# column 2: editing type
# column 3: sequence context. In the middle the edited base.
3430 AG CAG
3440 AG TAA
3441 AG AAT
3462 AG AAA
3456 AG TAC
3461 AG TAA
3481 AG CAG
</pre>
</div>
<li>SRR11550046.bwa.dis with the number of substitutions per each position along the read. It contains statistics for overlapping paired reads and the estimated error rate.</li>
<div class="highlight-python">
<pre>
#Example from overlap
>overlap
Reads in overlap: 4501245
Bases in overlap: 91212663
Errors in overlap: 17258
Error Rate: 0.000189
</pre>
</div> 
<li>SRR11550046.bwa.reads with the aligniment summary per read</li>
<div class="highlight-python">
<pre>
# Example from the content
# column 1: read ID
# column 2: first (1) or second (2) read of the pair
# columns 3 and 4: start and end coordinates of the alignment
# columns 5 and 6: start and end positions of the aligned portion
# column 7: percentage of aligned read
# column 8: comma separated list of genomic positions with variants
# column 9: comma separated list of read positions with variants
# column 10: comma separated list of variants
# column 11: number of variants with a phred score lower than the prefixed value (30)
# column 12: inferred strand
# column 13: 1: read is reverse 0: read is forward
# column 14: 1: mean phred score quality of variants
S100007198L1C004R0240315765 1 31 78 1 48 62.3 64,71 34,41 CT,CT 0 + 1 32.0
S100007198L1C005R0420577409 1 31 76 1 46 59.7 0 0 0 0 + 0
S100007198L1C001R0050400429 2 32 108 1 77 100.0 90 59 GA 0 + 0 32.0
S100007198L1C003R0060067932 1 29852 29903 1 52 67.5 29889,29891 38,40 AG,AC 0 + 1 40.0
S100007198L1C003R0580499073 1 29854 29903 2 51 64.9 29868,29869 15,16 GC,AT 0 + 1 31.5
</pre>
</div> 
<li>SRR11550046.bwa.cls with the list of reads with hyper editing in the same format of the .reads file</li>
<div class="highlight-python">
<pre>
# Example from the content
# See above for column explanation
S100007198L1C006R0100282147 2 11453 11529 1 77 100.0 11474,11476,11477,11486,11492,11493,11510,11513 22,24,25,34,40,41,58,61 AG,AG,AG,AG,AG,AG,AG,AG 0 + 0 38.375
S100007198L1C005R0580712041 2 12594 12670 1 77 100.0 12608,12611,12618,12650 15,18,25,57 AG,AG,AG,AG 0 + 0 39.25
S100007198L1C001R0020717137 1 12712 12788 1 77 100.0 12737,12740,12745,12777,12783 26,29,34,66,72 AG,AG,AG,AG,AG 0 + 1 37.6
S100007198L1C005R0580712041 1 12727 12803 1 77 100.0 12777,12788,12796,12800 51,62,70,74 AG,AG,AG,AG 0 + 1 40.0
S100007198L1C002R0450163396 1 16025 16101 1 77 100.0 16030,16032,16040,16061,16100 6,8,16,37,76 AG,AG,AG,AG,AG 0 + 1 39.6
S100007198L1C004R0090424767 2 17233 17309 1 77 100.0 17236,17269,17270,17276,17292,17300 4,37,38,44,60,68 AG,AG,AG,AG,AG,AG 0 + 0 38.8333333333
</pre>
</div> 
<li>SRR11550046.bwa.sum.cls with context values and a merging of overlapping hyper edited reads</li>
<div class="highlight-python">
<pre>
# Example from the content
# context one base up and down the edited bases
A_up 0.398004
C_up 0.166297
G_up 0.084257
T_up 0.351441
A_up 0.359202
C_up 0.186253
G_up 0.222838
T_up 0.231707
# column 1: flag (>) to extract only merged clusters
# column 2: variant type
# columns 3 and 4: start and end coordinates of the merged cluster
# column 5: strand of the merged cluster
# column 6: length of merged cluster
# column 7: number of variants in the merged cluster
# column 8: comma separated list of edits in the merged cluster
> AG 3205 3263 + 59 5 3205,3210,3222,3238,3263
> AG 3414 3507 + 94 12 3414,3417,3424,3430,3440,3441,3456,3461,3462,3481,3503,3507
> AG 4043 4107 + 65 6 4102,4107,4043,4049,4050,4079
</pre>
</div> 
</ul>
/opt/exp_soft/biomed/epicardi/python/bin/python /lustrehome/epicardi/home/covid/runREDItools.py Vero_SCV2-0.bam


