# SARS-COV-2 RNA editing
<b>A-to-I RNA editing in SARS-COV-2</b>

<p>RNA editing by adenosine deamination is a key molecular mechanism to modulate the innate immune response. During the course of a viral infection, it can have antiviral effects by destabilizing double stranded RNAs (dsRNAs) or proviral effects by suppressing the immune response.<br>
Here we provide the computational workflow to detect A-to-I RNA editing in the SARS-COV-2 genome. Our scripts enable the identification of hyper-edited reads as well as A-to-I editing at single nucleotide resolution.</p><br>
<b>Requirements</b>
<ul>
<li><b>Software</b></li>
<ul>
<li><a href="https://github.com/OpenGene/fastp">FASTP</a> for read trimming</li>
<li><a href="http://bio-bwa.sourceforge.net/">BWA</a> for read mapping</li>
<li><a href="http://research-pub.gene.com/gmap/">GSNAP</a> for read mapping</li>
<li><a href="https://lomereiter.github.io/sambamba/">sambamba</a> for filtering of SAM/BAM files</li>
<li><a href="http://www.htslib.org/">SAMtools</a> for manipulating BAM files</li>
<li><a href="https://github.com/BioinfoUNIBA/REDItools">REDItools v2</a> for RNA editing calling</li>
<li><a href="https://gsl.hudsonalpha.org/information/software/bam2fastq">bam2fastq</a> to estract reads from BAM files</li>
<li><a href="https://www.python.org/downloads/source/">python 2.x</a> to run scripts and the following mandatory modules:</li>
  <ul>
    <li><a href="https://pysam.readthedocs.io/en/latest/api.html">Pysam</a></li>
  </ul>
</ul>
<li><b>Genomic Data</b></li>
<ul>
<li>SARS-COV-2 genome (<a href="https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2">NC045512.2</a>)</li>
<li>Human Genome <a href="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz">hg19</a></li>
<li>Reads from <a href="https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA625518">PRJNA625518</a> BioProject - RNAseq from infected Calu-3 cells (only paired end runs):</li>
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
<p>The folder structure should be:</p>
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
<p><b>Detect hyper-edited reads</b></p>
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
<li>.bwa.Badreads with a list of low quality reads</li>
<li>SRR11550046.bwa.context with a list of nucleotides surrounding the editing sites to calculate the sequence context</li>
<div class="highlight-python">
<pre>
#Example of the content
# column 1: genomic position
# column 2: editing type
# column 3: sequence context.SRR11550046 In the middle the edited base.
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
<p><b>Profiling of RNA editing at single nucleotide level</b></p>
<p>Run REDItools v2</p>
<div class="highlight-python">
<pre>
reditools.py -r NC045512.fa -o SRR11550046.edi -f SRR11550046.bam -m myhomo -os 4 -q 30 -bq 30 -l 0 -s 2
</pre>
</div>
<p>Compare BWA and GSNAP alignments</p>
<div class="highlight-python">
<pre>
# run SubstitutionsPerSequence.py on gsnap BAM file
python SubstitutionsPerSequence.py SRR11550046.gsnap.bam 2 n
# compare aligned reads by CompareReads.py. It takes as input:
# BWA.reads file 
# GSNAP.reads file
# prefix for output files
python CompareReads.py SRR11550046.bwa.bam SRR11550046.gsnap.bam Comparison
</pre>
</div>
<p>CompareReads.py outputs 2 files:</p>
<ul>
<li>Comparison.badPos with a list of positions to exclude</li>
<li>Comparison.badReads with a list of reads to exclude</li>
</ul>
<p>Correct variants detected by REDItools</p>
<div class="highlight-python">
<pre>
# run corr.py. It takes as input:
# BAM file from bwa
# reditools output file to correct
# reference fasta file of viral genome
# strand info 1: FR Second Strand 2: FR First Strand 0: Unstranded
# file with positions to exclude (.badPos file)
# file with reads to exclude (.badReads file). Multiple files can be passed, separated by commas
# it requires also the list of known variants - knownVariants.txt - (modify its path inside the script)
python corr.py SRR11550046.bwa.bam SRR11550046.edi NC045512.fa 2 Comparison.badPos Comparison.badReads,SRR11550046.bwa.Badreads,SRR11550046.gsnap.Badreads
</pre>
</div>
<p>corr.py outputs 1 file:</p>
<ul>
<li>SRR11550046.corr.edi with a list of corrected positions using the REDItools output format</li>
<div class="highlight-python">
<pre>
# Example from the content
NC_045512.2	48	C	1	20465	37.71	[0, 20460, 0, 5]	CT	0.000244319569998	-	-	-	-	-
NC_045512.2	50	C	1	22437	37.60	[0, 22430, 0, 7]	CT	0.000311984668182	-	-	-	-	-
NC_045512.2	53	G	1	26391	39.10	[6, 0, 26383, 2]	GA	0.000227367463716	-	-	-	-	-
NC_045512.2	56	G	1	27791	39.07	[1, 0, 27786, 4]	GT	0.000143936667866	-	-	-	-	-
NC_045512.2	59	C	1	26491	37.73	[1, 26485, 1, 4]	CT	0.000151006077995	-	-	-	-	-
</pre>
</div> 
</ul>
<p>Select RNA editing positions using the selectPositions.py script from the REDItools package.</p>
<div class="highlight-python">
<pre>
# select positions using as minimal allele frequency two times the error rate from read overlap
# Error Rate: 0.000189 - Cut-off value: 0.000378
selectPositions.py -i SRR11550046.corr.edi -c 20 -v 4 -f 0.000378 -o SRR11550046.corr.sel.err.edi
# Calculate the distribution of detected variants
python subCount2.py SRR11550046.corr.sel.err.edi > SRR11550046.corr.sel.err.edi.distro
# Example of the content
# Row 1: variant types
# Row 2: number of detected variants per type
# Row 3: total number of detected variants
# Row 4: Occurrence of each variant type (in %)
AC	AG	AT	CA	CG	CT	GA	GC	GT	TA	TC	TG
8	62	15	5	2	235	14	14	23	21	143	22
564	564	564	564	564	564	564	564	564	564	564	564
1.41843971631	10.9929078014	2.65957446809	0.886524822695	0.354609929078	41.6666666667	2.48226950355	2.48226950355	4.0780141844	3.72340425532	25.3546099291	3.90070921986
</pre>
</div>
<div class="section" id="contact">
<b>Contact</b>
<ul class="simple">
<li>Ernesto Picardi: ernesto.picardi@gmail.com</li>
</ul>
</div>

<br><br><br><br><br>
Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg


This work is licensed under the Creative Commons Attribution - Non-Commercial 4.0 International License (CC BY-NC 4.0).
For any commercial use or license request for market activities, interested parties are invited to contact the Technology Transfer Office (TTO) of the University of Bari Aldo Moro, copyright holder. For a copy of the license, please visit https://creativecommons.org/licenses/by-nc/4.0/
