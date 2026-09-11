Lets get started with training a ChromBPNet model with a pre-trained bias model on the [downloaded](https://github.com/kundajelab/chrombpnet/wiki/Download-test-data) and [preprocessed](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing) tutorial data.

## Step 1

We will first start by downloading a pre-trained bias model provided with this github repo [here](https://zenodo.org/records/7443683/files/bias_models.zip?download=1). 

```
mkdir ~/chrombpnet_tutorial/bias_model
wget https://storage.googleapis.com/chrombpnet_data/input_files/bias_models/ATAC/ENCSR868FGK_bias_fold_0.h5 -O ~/chrombpnet_tutorial/bias_model/ENCSR868FGK_bias_fold_0.h5
```

- TODO: Add notes on how to pick a bias model

## Step 2

Use the pre-trained bias model to train a bias-factorized ChromBPNet model using the command below

```
chrombpnet pipeline \
        -ibam ~/chrombpnet_tutorial/data/downloads/merged.bam \
        -d "ATAC" \
        -g ~/chrombpnet_tutorial/data/downloads/hg38.fa \
        -c ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes \
        -p ~/chrombpnet_tutorial/data/peaks_no_blacklist.bed \
        -n ~/chrombpnet_tutorial/data/output_negatives.bed \
        -fl ~/chrombpnet_tutorial/data/splits/fold_0.json \
        -b ~/chrombpnet_tutorial/bias_model/ENCSR868FGK_bias_fold_0.h5 \
        -o ~/chrombpnet_tutorial/chrombpnet_model/
```
The command above outputs quality check report in two different formats - html and pdf. For your convenience and reference we provide links to both the outputs here - [html](https://storage.googleapis.com/chrombpnet_data/chrombpnet_model/evaluation/k562_overall_report.html) and [pdf](https://storage.googleapis.com/chrombpnet_data/chrombpnet_model/evaluation/k562_overall_report.pdf).

>**Important Notes:** 

> Every time you train the chrombpnet model please read the reports carefully to understand what is expected of the model and how you can correct if the bias correction is incomplete. 

> The two most important factors to consider while training a chrombpnet model are that the chrombpnet_nobias.h5 (1) learns only the expected TF motifs and (2) it does not learn any enzyme bias. This will be reflected in the TFModisco and marginal footprinting. 

> If the above criteria are not met, then try using a different bias model.

For general usage of this command you can run `chrombpnet pipeline -h` or refer to the documentation below. This command is intended to train and do quality checks on the chrombpnet model. You also have the option of performing the entire pipeline in two individual commands `chrombpnet train` and `chrombpnet qc`.

### Usage

```
chrombpnet pipeline [-h] -g GENOME -c CHROM_SIZES (-ibam INPUT_BAM_FILE | -ifrag INPUT_FRAGMENT_FILE | -itag INPUT_TAGALIGN_FILE) -o OUTPUT_DIR -d {ATAC,DNASE} -p PEAKS -n NONPEAKS -fl CHR_FOLD_PATH [-oth OUTLIER_THRESHOLD] [--ATAC-ref-path ATAC_REF_PATH] [--DNASE-ref-path DNASE_REF_PATH] [--num-samples NUM_SAMPLES] [-il INPUTLEN] [-ol OUTPUTLEN] [-s SEED] [-e EPOCHS] [-es EARLY_STOP] [-l LEARNING_RATE] [-track [TRACKABLES [TRACKABLES ...]]] [-a ARCHITECTURE_FROM_FILE] [-fp FILE_PREFIX] -b BIAS_MODEL_PATH [-sr NEGATIVE_SAMPLING_RATIO] [-fil FILTERS] [-dil N_DILATION_LAYERS]  [-j MAX_JITTER] [-bs BATCH_SIZE]
```

### Input Format

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
                        10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) +
                        summit (10th col).
  -n NONPEAKS, --nonpeaks NONPEAKS
                        10 column bed file of non-peak regions, centered at summit (10th column)
  -fl CHR_FOLD_PATH, --chr-fold-path CHR_FOLD_PATH
                        Fold information - dictionary with test,valid and train keys and values with corresponding
                        chromosomes
  -b BIAS_MODEL_PATH, --bias-model-path BIAS_MODEL_PATH
                        Path for a pretrained bias model

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
                        list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,
                        val_logcount_predictions_loss,val_loss,val_profile_predictions_loss
  -a ARCHITECTURE_FROM_FILE, --architecture-from-file ARCHITECTURE_FROM_FILE
                        Model to use for training
  -fp FILE_PREFIX, --file-prefix FILE_PREFIX
                        File prefix for output to use. All the files will be prefixed with this string if provided.
  -sr NEGATIVE_SAMPLING_RATIO, --negative-sampling-ratio NEGATIVE_SAMPLING_RATIO
                        Ratio of negatives to positive samples per epoch
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

#### Output Format

The output directory will be populated as follows -

```
models\
	bias_model_scaled.h5
	chrombpnet.h5
	chrombpnet_nobias.h5 (TF-Model i.e model to predict bias corrected accessibility profile) 
logs\
	chrombpnet.log (loss per epoch)
	chrombpnet.log.batch (loss per batch per epoch)
	(..other hyperparameters used in training)
	
auxilary\
	filtered.peaks
	filtered.nonpeaks
	...

evaluation\
	overall_report.pdf
	overall_report.html
	bw_shift_qc.png 
	bias_metrics.json 
	chrombpnet_metrics.json
	chrombpnet_only_peaks.counts_pearsonr.png
	chrombpnet_only_peaks.profile_jsd.png
	chrombpnet_nobias_profile_motifs.pdf
	chrombpnet_nobias_counts_motifs.pdf
	chrombpnet_nobias_max_bias_response.txt
	chrombpnet_nobias.....footprint.png
	...
```


- `overall_report.html` is an html summary of both the training and motifs learnt by the chrombpnet model. It also has guidance on quality checks to consider for the chrombpnet model. This report references all the images in the `evaluation/` folder.
- `overall_report.pdf` is the pdf rendered file of the html report.
- For a full description of the remaining files and folders refer to the [next section](https://github.com/kundajelab/chrombpnet/wiki/Output-format) on output formats. 
- If the bias model transfer has failed try a different bias model or train a custom bias model following the instructions in [here](https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training).
- The motif list in `overall_report.pdf` and `overall_report.html` is not representative of the motif list on all the peaks. This motif list is obtained from a random subset of peaks (30,000) for bias quick quality check. In order to get the representative motif summary from TF-Modisco on all enhancers please fetch [contribution scores](https://github.com/kundajelab/chrombpnet/wiki/Generate-contribution-score-bigwigs) on all peaks and utilize [TFModisco lite](https://github.com/kundajelab/chrombpnet/wiki/Denovo-motif-discovery) with `-n` set to 1 million. 


