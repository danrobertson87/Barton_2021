### Need to set
# PATH_TO/Reference.chrom.sizes
##

mkdir calibrate

# calculate OR value
# Get the number of mapped reads from flagstat

# calibration input
num1=`cat  miniMap2_out_calibration_genome_only/$1.calibration.only.mm2.chrMT.bam.fs | grep N/A|grep mapped | grep N/A |cut -f1 -d ' '`
echo "$1 pombe: $num1 (Wc)"

# experimental IP
num2=`cat  miniMap2_out_experimental_genome_only/$2.experimental.only.mm2.MT.rDNA.bam.fs | grep N/A|grep mapped | grep N/A |cut -f1 -d ' '`
echo "$2 sacCer3: $num2 (IPx)"

# calibration IP
num3=`cat  miniMap2_out_calibration_genome_only/$2.calibration.only.mm2.chrMT.bam.fs | grep N/A|grep mapped | grep N/A |cut -f1 -d ' '`
echo "$2 pombe: $num3 (IPc)"

# experimental input
num4=`cat  miniMap2_out_experimental_genome_only/$1.experimental.only.mm2.MT.rDNA.bam.fs | grep N/A|grep mapped | grep N/A |cut -f1 -d ' '`
echo "$1 sacCer3: $num4 (Wx)"

WcIPx=`perl -e "print $num1*$num2"`
WxIPc=`perl -e "print $num4*$num3"`

OR=`perl -e "print $WcIPx / $WxIPc;"`
echo "The occupancy ratio is $OR"

# Create BigWigs (Calibrated)
cat miniMap2_out_experimental_genome_only/$2.experimental.only.MT.rDNA.rpm.bedgraph | awk '{print $1 "\t" $2 "\t" $3 "\t" $4*'$OR'}' > calibrate/$2.experimental.only.MT.rDNA.rpm.calibrated.bedgraph

# Convert to bigwig
wigToBigWig calibrate/$2.experimental.only.MT.rDNA.rpm.calibrated.bedgraph PATH_TO/Reference.chrom.sizes calibrate/$2.experimental.only.MT.rDNA.rpm.calibrated.bw