Here we discuss different preprocessing steps using the example source data previously [downloaded](https://github.com/kundajelab/chrombpnet/wiki/Download-test-data). 

## Mapping

Merge isogenic replicates of experiments together to yield consolidated reads. To merge bam files we use `samtools merge`. 

```
# merge and index bam files

samtools merge -f ~/chrombpnet_tutorial/data/downloads/merged_unsorted.bam ~/chrombpnet_tutorial/data/downloads/rep1.bam  ~/chrombpnet_tutorial/data/downloads/rep2.bam  ~/chrombpnet_tutorial/data/downloads/rep3.bam
samtools sort -@4 ~/chrombpnet_tutorial/data/downloads/merged_unsorted.bam -o ~/chrombpnet_tutorial/data/downloads/merged.bam
samtools index ~/chrombpnet_tutorial/data/downloads/merged.bam
```

- Fragment files obtained from Chromap, Cellranger and similar tools are suitable with ChromBPNets pipeline. For single-cell analyses, you may wish to aggregate fragments by cluster instead of a sample to train cluster-specific ChromBPNet models.
- TagAlign files are also suitable with ChromBPNets pipeline.
- When using the ChromBPNet pipelines to train models, you do not need to worry about pre-shifting the fragments, bams or tagAligns. The training tools are designed to automatically detect and correct any shifts in the files, specifically for ATAC and DNase, so you can run training without any additional shifting.

## Filtering

- If your bams are of mixed types (i.e includes both paired-end and single-end replicates) do filtering before the above step (merging). 
- Filter the input files for quality metrics  and remove duplicates appropriately. When working with paired-end data, duplicates can be easily identified and removed because we have paired-end reads. However, when working with single-end reads, it is not possible to accurately identify duplicates, so it is not recommended to remove duplicates for single-end data as it might result in throwing away actual signal information.
- The input files [ENCSR868FGK](https://github.com/kundajelab/chrombpnet/wiki/Download-test-data) for ATAC-seq paired-end data used in the tutorial have already been appropriately filtered and duplicates have been handled using the ENCODE ATAC-seq protocol, so no further action is required for the tutorial. Therefore, you don't have to do anything with the input files for the purpose of this tutorial.
- If you are downloading the data from [ENCODE portal](https://www.encodeproject.org/), download the filtered bams for paired-end data and  unfiltered bams for single-end data. For filtered bams no further filtering or de-duplication is required. For the unfiltered bams filter for quality metrics using (`samtools view -b -@50 -F780 -q30`).
 
## Peak-calling

- We recommend using relaxed peak-calls (MACS2 with p-value parameter threshold set to 0.01) similar to the overlap peaks specifications in [ENCODE ATAC-seq protocol](https://github.com/ENCODE-DCC/atac-seq-pipeline).
- For the tutorial experiment ENCSR868FGK relaxed peak call are already provided on ENCODE portal and can  be downloaded as follows - 

```
# download overlap peaks (default peaks on ENCODE)
wget https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz -O ~/chrombpnet_tutorial/data/downloads/overlap.bed.gz
```

- If you are downloading the data from the [ENCODE portal](https://www.encodeproject.org/) you can download the peaks flagged default for ATAC-seq datasets. For DNASE-seq datasets you might have to use MACS2 and follow the ENCODE ATAC-seq protocol for peak-calling.


- Ensure that the peak regions do not intersect with the blacklist regions (extended by 1057 bp on both sides) by running the command below. 

```
bedtools slop -i ~/chrombpnet_tutorial/data/downloads/blacklist.bed.gz -g ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes -b 1057 > ~/chrombpnet_tutorial/data/downloads/temp.bed
bedtools intersect -v -a ~/chrombpnet_tutorial/data/downloads/overlap.bed.gz -b ~/chrombpnet_tutorial/data/downloads/temp.bed  > ~/chrombpnet_tutorial/data/peaks_no_blacklist.bed
```

## Define train, validation and test chromosome splits

To train ChromBPNet models, we divide the dataset into three groups by chromosomes: training, validation, and testing. We use standard ENCODE cross-validation folds for human chromosomes in our work. This tutorial demonstrates how to generate the json file for fold 0 (`~/chrombpnet_tutorial/data/splits/fold_0.json`) and train the model on that fold. The remaining folds for human chromosomes can be found here - ([fold_0](https://zenodo.org/record/7445373/files/fold_0.json?download=1), [fold_1](https://zenodo.org/record/7445373/files/fold_1.json?download=1), [fold_2](https://zenodo.org/record/7445373/files/fold_2.json?download=1), [fold_3](https://zenodo.org/record/7445373/files/fold_3.json?download=1), [fold_4](https://zenodo.org/record/7445373/files/fold_4.json?download=1)).

>**Note:** 
> To create the folds we want to specify only the  training, validation and testing chromosomes. For the human chromosomes we use chr1 to chr22, chrX and chrY and filter the remaining.
```
head -n 24  ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes >  ~/chrombpnet_tutorial/data/downloads/hg38.chrom.subset.sizes
```

```
mkdir ~/chrombpnet_tutorial/data/splits
chrombpnet prep splits -c ~/chrombpnet_tutorial/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ~/chrombpnet_tutorial/data/splits/fold_0
```

For general usage of this command you can run `chrombpnet prep splits -h` or refer to the documentation below. Make sure that you only have chromosomes of interest (i.e chromosomes you want to include for training, validation and testing) in the `-c CHROM_SIZES` file.

### Usage 

```
chrombpnet prep splits [-h] -op OUTPUT_PREFIX -c CHROM_SIZES -tcr [TEST_CHROMS [TEST_CHROMS ...]] -vcr
                              [VALID_CHROMS [VALID_CHROMS ...]]
```

### Input Format

```
required arguments:
  -op OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Path prefix to store the fold information (appended with .json)
  -c CHROM_SIZES, --chrom-sizes CHROM_SIZES
                        TSV file with chromosome sizes. All chromosomes from the first column of chrom sizes file are used
  -tcr [TEST_CHROMS [TEST_CHROMS ...]], --test-chroms [TEST_CHROMS [TEST_CHROMS ...]]
                        Chromosomes to use for test
  -vcr [VALID_CHROMS [VALID_CHROMS ...]], --valid-chroms [VALID_CHROMS [VALID_CHROMS ...]]
                        Chromosomes to use for validation
```

### Output Format

- `output_prefix`.json: The script outputs a json file with the prefix `output_prefix`. 

>**Note:** 
>Note that prefix can include a directory path and prefix for the output file. Make sure that the directory in output_prefix exists.

## Generate non-peaks (background regions)

Here we will generate non-peak background regions (`~/chrombpnet_tutorial/data/output_negatives.bed`) that GC-match with the peak regions. We will use these regions in [ChromBPNet model training](https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training), [bias model training](https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training) and as background regions to get [marginal footprints](https://github.com/kundajelab/chrombpnet/wiki/Marginal-footprinting) in the upcoming sections of this tutorial. 

```
chrombpnet prep nonpeaks -g ~/chrombpnet_tutorial/data/downloads/hg38.fa -p ~/chrombpnet_tutorial/data/peaks_no_blacklist.bed -c  ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes -fl ~/chrombpnet_tutorial/data/splits/fold_0.json -br ~/chrombpnet_tutorial/data/downloads/blacklist.bed.gz -o ~/chrombpnet_tutorial/data/output
```

For general usage of this command you can run `chrombpnet prep nonpeaks -h` or refer to the documentation below. 

### Usage

```
chrombpnet prep nonpeaks [-h] -g GENOME -o OUTPUT_PREFIX -p PEAKS -c CHROM_SIZES -fl CHR_FOLD_PATH [-il INPUTLEN] [-st STRIDE]
                                [-npr NEG_TO_POS_RATIO_TRAIN] [-br BLACKLIST_REGIONS]
```

### Input Format

```
required arguments:
  -g GENOME, --genome GENOME
                        reference genome file
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        output BED file prefix to store the gc content of binned genome. suffix .bed will be appended by the code. If
                        the prefix contains a directory path make sure it exists.
  -p PEAKS, --peaks PEAKS
                        10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th
                        col).
  -c CHROM_SIZES, --chrom-sizes CHROM_SIZES
                        Chrom sizes file
  -fl CHR_FOLD_PATH, --chr-fold-path CHR_FOLD_PATH
                        Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes

optional arguments:
  -il INPUTLEN, --inputlen INPUTLEN
                        inputlen to use to make bins and find gc content. Default set to 2114. 
  -st STRIDE, --stride STRIDE
                        stride to use for shifting the bins. Default set to 1000.
  -npr NEG_TO_POS_RATIO_TRAIN, --neg-to-pos-ratio-train NEG_TO_POS_RATIO_TRAIN
                        Ratio of negatives to positives to sample in training set (test set always has 1:1 positive to negatives ratio). 
                        Default set to 2. 
  -br BLACKLIST_REGIONS, --blacklist-regions BLACKLIST_REGIONS
                        TSV file with 3 columns - chr, start, end. Backlisted regions not to be used for training.
  -s SEED
```

### Output Format

- `output_prefix`_negatives.bed: A bed format file with 10 columns that contains non-peak regions. These regions are selected to have a neg_to_pos_ratio_train ratio with the number of regions in peaks in the training and validation chromosomes and a 1:1 ratio with the number of regions in peaks in the test chromosomes.
- `output_prefix`_auxiliary: A directory containing intermediate files generated as a part of this script. Will be removed in future releases if considered redundant. 


>**Note:** 
> Note that prefix can include a directory path and prefix for the output file. Make sure that the directory in output_prefix exists.
> This script typically runs in 10 minutes for the tutorial data given here. For any other custom dataset and reference genome, if this script does not converge to an output in a long time (>20 minutes). Try reducing the optional argument`-st` to 500 or 250. In future release we will try to provide a better warning or an iterative solution for this. 


