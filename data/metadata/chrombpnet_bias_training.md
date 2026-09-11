To achieve good bias correction it is important to have a good and complete bias model. While bias models can be transferred across experiments done under similar settings (access a list of pretrained bias models [here](https://zenodo.org/records/7443683)), one might need to train a new bias model when the experiment conditions change. In this section of the tutorial we detail the command by which one can train their own bias model using again the [downloaded](https://github.com/kundajelab/chrombpnet/wiki/Download-test-data) and [preprocessed](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing) tutorial data.

```
chrombpnet bias pipeline \
        -ibam ~/chrombpnet_tutorial/data/downloads/merged.bam \
        -d "ATAC" \
        -g ~/chrombpnet_tutorial/data/downloads/hg38.fa \
        -c ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes \
        -p ~/chrombpnet_tutorial/data/peaks_no_blacklist.bed \
        -n ~/chrombpnet_tutorial/data/output_negatives.bed \
        -fl ~/chrombpnet_tutorial/data/splits/fold_0.json \
        -b 0.5 \
        -o ~/chrombpnet_tutorial/bias_model/ \
        -fp k562 \
```

The command above outputs quality check report in two different formats - html and pdf. For your convenience and reference we provide links to both the outputs here - [html](https://storage.googleapis.com/chrombpnet_data/bias_model/evaluation/k562_overall_report.html) and [pdf](https://storage.googleapis.com/chrombpnet_data/bias_model/evaluation/k562_overall_report.pdf). 

>**Important Notes:** 

> Every time you train the bias model please read the reports carefully to understand what is expected of the model and how you can correct if the bias model is inaccurate. 

> The three most important factors to consider while training a bias model are (1) it learns the expected bias motif (2) it does not learn any Transcription Factor motifs. If this condition is not met the bias model will regress out TF activity in addition to bias activity from accessibility profiles and (3) the bias model trained in non-peaks has similar GC bias as peaks. If this condition is not met the bias model transfer fails when used in peaks.

>The bias threshold factor is a single hyperparameter that will help you control all the above three factors. For starters we recommend using 0.5 for "ATAC" assay type and 0.8 for "DNASE" assay type. To understand in detail how to tune this hyperparameter refer to [FAQ](https://github.com/kundajelab/chrombpnet/wiki/FAQ) or Input Format section below. 

For general usage of this command you can run `chrombpnet bias pipeline -h` or refer to the documentation below. This command is intended to train and do quality checks on the bias model. You also have the option of performing the entire pipeline in two individual commands `chrombpnet bias train` and `chrombpnet bias qc`. 
 

### Usage 

```
chrombpnet bias pipeline [-h] -g GENOME -c CHROM_SIZES (-ibam INPUT_BAM_FILE | -ifrag INPUT_FRAGMENT_FILE | -itag INPUT_TAGALIGN_FILE) -o OUTPUT_DIR -d {ATAC,DNASE} -p PEAKS -n NONPEAKS -fl CHR_FOLD_PATH [-oth OUTLIER_THRESHOLD] [--ATAC-ref-path ATAC_REF_PATH] [--DNASE-ref-path DNASE_REF_PATH] [--num-samples NUM_SAMPLES] [-il INPUTLEN] [-ol OUTPUTLEN] [-s SEED] [-e EPOCHS] [-es EARLY_STOP] [-l LEARNING_RATE] [-track [TRACKABLES [TRACKABLES ...]]] [-a ARCHITECTURE_FROM_FILE] [-fp FILE_PREFIX] -b BIAS_THRESHOLD_FACTOR [-fil FILTERS] [-dil N_DILATION_LAYERS] [-j MAX_JITTER] [-bs BATCH_SIZE]
```

#### Input Format

```
required arguments:
  -g GENOME, --genome GENOME
                        reference genome fasta file
  -c CHROM_SIZES, --chrom-sizes CHROM_SIZES
                        Chrom sizes file
  -ibam INPUT_BAM_FILE, --input-bam-file INPUT_BAM_FILE
                        Input BAM file
  -ifrag INPUT_FRAGMENT_FILE, --input-fragment-file INPUT_FRAGMENT_FILE
                        Input fragment file
  -itag INPUT_TAGALIGN_FILE, --input-tagalign-file INPUT_TAGALIGN_FILE
                        Input tagAlign file
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output dir (path/to/output/dir)
  -d {ATAC,DNASE}, --data-type {ATAC,DNASE}
                        assay type
  -p PEAKS, --peaks PEAKS
                        10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).
  -n NONPEAKS, --nonpeaks NONPEAKS
                        10 column bed file of non-peak regions, centered at summit (10th column)
  -fl CHR_FOLD_PATH, --chr-fold-path CHR_FOLD_PATH
                        Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes
  -b BIAS_THRESHOLD_FACTOR, --bias-threshold-factor BIAS_THRESHOLD_FACTOR
                        A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x 
                        min(count over peakregions). Recommended start value 0.5 for ATAC and 0.8 for DNase.

optional arguments:
  -oth OUTLIER_THRESHOLD, --outlier-threshold OUTLIER_THRESHOLD
                        threshold to use to filter outlies
  --ATAC-ref-path ATAC_REF_PATH
                        Path to ATAC reference motifs (ATAC.ref.motifs.txt used by default)
  --DNASE-ref-path DNASE_REF_PATH
                        Path to DNASE reference motifs (DNASE.ref.motifs.txt used by default)
  --num-samples NUM_SAMPLES
                        Number of reads to sample from BAM/fragment/tagAlign file for shift estimation
  -il INPUTLEN, --inputlen INPUTLEN
                        Sequence input length
  -ol OUTPUTLEN, --outputlen OUTPUTLEN
                        Prediction output length
  -s SEED, --seed SEED  seed to use for model training
  -e EPOCHS, --epochs EPOCHS
                        Maximum epochs to train
  -es EARLY_STOP, --early-stop EARLY_STOP
                        Early stop limit, corresponds to 'patience' in callback
  -l LEARNING_RATE, --learning-rate LEARNING_RATE
                        Learning rate for model training
  -track [TRACKABLES [TRACKABLES ...]], --trackables [TRACKABLES [TRACKABLES ...]]
                        list of things to track per batch, such as
                        logcount_predictions_loss, loss,profile_predictions_loss, val_logcount_predictions_loss, val_loss, 
                        val_profile_predictions_loss
  -a ARCHITECTURE_FROM_FILE, --architecture-from-file ARCHITECTURE_FROM_FILE
                        Model to use for training
  -fp FILE_PREFIX, --file-prefix FILE_PREFIX
                        File prefix for output to use. All the files will be prefixed with this string if provided.
  -fil FILTERS, --filters FILTERS
                        Number of filters to use in chrombpnet mode
  -dil N_DILATION_LAYERS, --n-dilation-layers N_DILATION_LAYERS
                        Number of dilation layers to use in chrombpnet model
  -j MAX_JITTER, --max-jitter MAX_JITTER
                        Maximum jitter applied on either side of region (default 500 for chrombpnet model)
  -bs BATCH_SIZE, --batch-size BATCH_SIZE
                        batch size to use for model training
```

- Only one of `-ibam`, `-ifrag` and `itag` can be used as arguments. Example files for supported types are provided here for reference - [bam](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_merged.bam), [fragment](https://storage.googleapis.com/chrombpnet_data/input_files/example.fragments.tsv), [tagalign](https://storage.googleapis.com/chrombpnet_data/input_files/example.tagAlign) 
- The recommended starting value to use for the `-b` argument or the `bias threshold factor is 0.5 for "ATAC" assay type and 0.8 for "DNASE" assay type. 


#### Output Format

The output directory will be populated as follows -

```
models\
	bias.h5
logs\
	bias.log (loss per epoch)
	bias.log.batch (loss per batch per epoch)
	(..other hyperparameters used in training)
	
intermediates\
	...

evaluation\
        overall_report.html
        overall_report.pdf
	pwm_from_input.png
        k562_epoch_loss.png 
	bias_metrics.json
	bias_only_peaks.counts_pearsonr.png
	bias_only_peaks.profile_jsd.png
	bias_only_nonpeaks.counts_pearsonr.png
	bias_only_nonpeaks.profile_jsd.png
        bias_predictions.h5
	bias_profile.pdf
	bias_counts.pdf
	...
```

- `overall_report.html` is an html summary of both the training and motifs learnt by the bias model. It also has guidance on quality checks to consider for the bias model. This report references all the remaining images in the `evaluation/` folder.
- `overall_report.pdf` is the pdf rendered file of the html report.
- For a full description of the remaining files and folders refer to the [next section](https://github.com/kundajelab/chrombpnet/wiki/Output-format) on output formats. 
