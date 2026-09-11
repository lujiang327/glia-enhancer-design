### I am having trouble loading chrombpnet.h5 file, what should I do?

If you are having trouble loading `chrombpnet.h5` model and you are seeing the following error `ValueError: bad marshal data (unknown type code)` its likely because the tensorflow version you are using is different from the tensorflow version used to train the model.

if you have `saved_model.pb` within a `chrombpnet/` directory I would recommend you load the model as follows `load_model('chrombpnet/')`. This model format should be portable across most tensorflow versions. 

If the error persists, run the following script [here](https://github.com/kundajelab/chrombpnet/tree/master/chrombpnet/helpers/postprocessing) on the files `chrombpnet_nobias.h5` and `bias_model_scaled.h5`. This script will recompile the `chrombpnet.h5` model for that tensorflow version you are using!

### How can I make prediction and contribution score bigwigs if I have a bed3 format? How can I convert it to NarrowPeak format to make it compatible with the code?

(1) Merge the input bed file (2) Make overlapping windows for the merged bed file using bedtools makewindows. (window size 1000bp and overlapping by at least 250bp so step size is 750 bp) (3) Create summits centered at the windows to make predictions. So 10th col will be (end-start)/2 and the remaining cols will be `.`

### All of the functions in the tutorial use fold_0, however, in general practice would you build multiple models across multiple folds?

Yes you will need to train a model for each fold. And then average the output bigwigs (both predictions and contribution scores) across folds. You will also need to average the contribution score h5's before running TFModisco. We will eventually release a simple tool for this in future. But for now you can average the bigwigs using existing tools (e.g. `wiggletools mean`) or average the `h5s` by loading them in python.


### What CUDA packages are relevant for this repository? 

This repository currently uses tensorflow 2.8.0. Documented in the [requirements.txt](https://github.com/kundajelab/chrombpnet/blob/master/requirements.txt). To find the relevant CuDNN and CUDA versions refer to the documentation [here](https://www.tensorflow.org/install/source#tested_build_configurations). For tensorflow 2.8.0 its `cudnn/8.1` and `cuda/11.2.0`. 
 
### (1) How to pick a pre-trained bias model? 

Pick a pre-trained bias model that is similar in experimental setting to the current dataset. For example, if your current dataset follows a Omni ATAC-seq protocol use a bias model trained on any Omni ATAC-seq protocol dataset. If multiple Omni ATAC-seq pre-trained bias model exists pick a bias model trained on the the closest biosample or highest read depth. 

### (2) What is the intuition in choosing the hyperparameter for `bias_threshold_factor`? How to retrain the bias model based on this ? 

Non peak regions used in bias model training are filtered based on The  `bias_threshold_factor` which is used as follows. The regions with total counts greater than 0.01_quantile(total counts in peaks)*bias_threshold_factor are filtered out. 

- If `bias_threshold_factor` is set to very low value you might filter out a lot of non-peak regions and the resulting training set might have GC distribution that is different from that in peak regions. This will result in bias model transfer because the model will learn a AT rich bias.

- If `bias_threshold_factor` is set to very high you might include non-peak regions with high counts. High counts are typically caused by Trnascription Factor (TF) motifs and this might lead to the bias model capturing cell-type specific TFs (which is not ideal as we want to regress out only the bias motifs effect and not cell-type specific TF motifs effect).

In choosing a value for `bias_threshold_factor` we recommend starting from 0.5 for ATAC and 0.8 for DNase. If these models capture TF motifs reduce the `bias_threshold_factor` by 0.1 and keep doing this until you can train a model that captures only bias motifs (without any TF motifs). If these models have high negative pearsonr (> -0.3) in peak regions but capture no TF motifs, increase the `bias_threshold_factor` by 0.1 until you  you can increase the pearsonr value (to a value > -0.3) while also not capturing any TF motifs. 

### (3) Why does the quality check reports involve getting contribution scores and TFModisco motifs only on subset of peaks?

It is computationally expensive to run these steps on all the peaks. Since we are looking for a quick sanity check for our models we will subsample 30K peaks from our original peak set for interpretation and TFMOdisco. The users are encouraged to eventually run `chrombpet contribs_bw` and `chrombpnet modisco_motifs` step on all the peak regions for interpretation and motif discovery.


### (4) How can I run ChromBPNet on Stanford sherlock cluster? 

[Sherlock](https://www.sherlock.stanford.edu/) does not support docker, but it does support [singularity](https://www.sherlock.stanford.edu/docs/software/using/singularity/). Singularity is installed by default on all sherlock nodes.
To run chrombpnet on sherlock: 

1) Make sure you are on a gpu node. The command below will start an interactive node session that will remain active for a day. For more information see documentation on [interactive nodes](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/#compute-nodes) and [srun](https://slurm.schedmd.com/srun.html). For doing scalable runs on sherlock we recommend using [sbatch](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/#batch-jobs) over srun.
```
srun -p gpu -c 4 --gres gpu:1 --time 1-0 --pty bash
```

2) Load appropriate cuDNN and CUDA versions when using srun and local condo setup. For example - for tensorflow 2.8.0 you will need to do `ml cudnn/8.1` and `ml cuda/11.2.0`. 

3) Use the `singularity exec` command, binding your data directory to a directory within the singularity container. An example of exec `chrombpnet pipeline` for the [tutorial example](https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training) is provided below: 

```
# we will mount all the required files in our current directory (`~/chrombpnet_tutorial/`), which contains the data and output folders from tutorial, to the `/mnt` directory on the singularity container and run chrombpnet pipeline by using the command below.

singularity  exec --nv -e --no-mount hostfs --bind ~/chrombpnet_tutorial/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet pipeline -ibam /mnt/data/downloads/merged.bam -d ATAC -g /mnt/data/downloads/hg38.fa -c /mnt/data/downloads/hg38.chrom.sizes -p /mnt/data/peaks_no_blacklist.bed -n /mnt/data/negatives.bed -f /mnt/data/splits/fold_0.json -b /mnt/bias_model/ENCSR868FGK_bias_fold_0.h5 -o /mnt/chrombpnet_model/
```

Note that the flags used with singularity exec are important. This is what each flag means: 

* `--nv` --> run on a GPU 


* `-e` --> clean environment before running the container 

* `--no-mount hostfs` --> do not mount the host filesystem in the singularity container. This avoids conflicting installations of python, and overwriting directories in the singularity image (i.e. avoids mounting `/scratch` on sherlock to `/scratch` in the singularity container, which we don't want as that is where the chrombpnet source files live in the container) 

* `--bind` --> binds paths on the host machine to paths in the container. The syntax is path1_on_host:path1_on_container,path2_on_host:path2_on_container,pathN_on_host:pathN_on_container. 