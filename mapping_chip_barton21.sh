### Need to set
#PATH/calibration_genome_reference.fa
#PATH/experimental_genome_reference.fa
#PATH/calibration_genome_reference.chrom.sizes
#PATH/experimental_genome_reference.chrom.sizes
###

### First mapping to **Calibration** Genome

mkdir miniMap2_out_calibration_genome

# Map raw/trimmed reads to calibration genome with Minimap2 (paired end)
minimap2 -ax sr PATH/calibration_genome_reference.fa cutadapt/$1.trimmed.fastq.gz cutadapt/$2.trimmed.fastq.gz > miniMap2_out_calibration_genome/$3.calibration.mm2.sam

# SAM to BAM to sort (unmapped reads only)
samtools view -Sbh -f 4 miniMap2_out_calibration_genome/$3.calibration.mm2.sam -o miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.bam

# sort unmapped calibration BAM (-n by qname)
samtools sort -n -@ 5 miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.bam -o miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.qsort.bam

# interleaved (keep in single fastq)
samtools bam2fq miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.qsort.bam > miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.fastq
gzip miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.fastq

# clean up 
rm miniMap2_out_calibration_genome/$3.calibration.mm2.sam
rm miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.bam 
rm miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.qsort.bam


### First mapping to **Experimental** Genome

mkdir miniMap2_out_experimental_genome

# Map raw/trimmed reads to experimental genome with Minimap2 (paired end)
minimap2 -ax sr PATH/experimental_genome_reference.fa cutadapt/$1.trimmed.fastq.gz cutadapt/$2.trimmed.fastq.gz > miniMap2_out_experimental_genome/$3.experimental.mm2.sam

# SAM to BAM to sort (unmapped reads only)
samtools view -Sbh -f 4 miniMap2_out_experimental_genome/$3.experimental.mm2.sam -o miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.bam

# sort unmapped sacer BAM (-n by qname)
samtools sort -n -@ 5 miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.bam -o miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.qsort.bam

# interleaved (keep in single fastq)
samtools bam2fq miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.qsort.bam > miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.fastq
gzip miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.fastq

# clean up
rm miniMap2_out_experimental_genome/$3.experimental.mm2.sam
rm miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.bam
rm miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.qsort.bam 


### Second mapping to **Calibration** Genome

mkdir miniMap2_out_calibration_genome_only

# Map unmapped experimental genome reads to calibration genome with Minimap2 (interleaved)
# input miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.fastq.gz
minimap2 -ax sr PATH/calibration_genome_reference.fa miniMap2_out_experimental_genome/$3.experimental.mm2.unmapped.fastq.gz > miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.sam

# SAM to BAM to sort (mapped reads only)
samtools view -Sbh -F 4 miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.sam -o miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam
samtools sort -@ 5 miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam -o miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.sorted.bam

rm miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam
rm miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.sam
mv miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.sorted.bam miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam

samtools index miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam

# remove chrMT/MTR/AB
samtools view -b miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.bam I II III > miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam
samtools index miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam
samtools flagstat miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam > miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam.fs

###
# Create BigWigs (READS PER MILLION)
# Get the number of mapped reads from flagstat
num1=`cat  miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam.fs | grep N/A | grep mapped | grep N/A |cut -f1 -d ' '`
echo "$3: $num1"

# Convert to RPM
RPM=`perl -e "print 1000000/$num1"`
echo "The scaling factor is $RPM"

# Create a bedgraph file scaled to RPM 
genomeCoverageBed -ibam  miniMap2_out_calibration_genome_only/$3.calibration.only.mm2.chrMT.bam -bg -scale $RPM > miniMap2_out_calibration_genome_only/$3.calibration.only.chrMT.rpm.bedgraph

# Convert to bigwig 
wigToBigWig miniMap2_out_calibration_genome_only/$3.calibration.only.chrMT.rpm.bedgraph PATH/calibration_genome_reference.chrom.sizes miniMap2_out_calibration_genome_only/$3.calibration.only.chrMT.rpm.bw


### Second mapping to **Experimental** Genome

mkdir miniMap2_out_experimental_genome_only

# Map unmapped calibration genome to experimental genome with Minimap2 (interleaved)
# input miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.fastq
minimap2 -ax sr ~genomes/s.cerevisiae/sacCer3/sacCer3.fa miniMap2_out_calibration_genome/$3.calibration.mm2.unmapped.fastq.gz > miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.sam

# SAM to BAM to sort (mapped reads only)
samtools view -Sbh -F 260 -f 3 miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.sam -o miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam
samtools sort -@ 5 miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam -o miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.sorted.bam

rm miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam
rm miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.sam
mv miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.sorted.bam miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam

samtools index miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam

# exclude chrMT
samtools view -b miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.bam chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI > miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.bam

samtools index miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.bam
samtools flagstat miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.bam > miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.bam.fs

# Remove rDNA regions (chrXII saturated with reads)
bedtools intersect -v -abam miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.bam -b rDNA.bed > miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.rDNA.bam
samtools index miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.rDNA.bam

###
# Create BigWigs (READS PER MILLION)
# Get the number of mapped reads from flagstat
num1=`cat  miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.rDNA.bam.fs | grep N/A|grep mapped | grep N/A |cut -f1 -d ' '`
echo "$3: $num1"

# Convert to RPM
RPM=`perl -e "print 1000000/$num1"`
echo "The scaling factor is $RPM"

# Create a bedgraph file scaled to RPM
genomeCoverageBed -ibam  miniMap2_out_experimental_genome_only/$3.experimental.only.mm2.MT.rDNA.bam -bg -scale $RPM  > miniMap2_out_experimental_genome_only/$3.experimental.only.MT.rDNA.rpm.bedgraph

# Convert to bigwig
wigToBigWig  miniMap2_out_experimental_genome_only/$3.experimental.only.MT.rDNA.rpm.bedgraph PATH/experimental_genome_reference.chrom.sizes miniMap2_out_experimental_genome_only/$3.experimental.only.MT.rDNA.rpm.bw