## Reconstructing phylogenetic trees from genome-wide somatic mutations in clonal samples

The R script in this repository takes in matrices of variant depth and total depth across from multiple clonal samples of the same donor, filters out germline and artefactual variants, constructs the phylogenetic tree topology and maps somatic mutations to the tree branches using the [TreeMut](https://github.com/nangalialab/treemut) package. 

For any queries, raise an issue on the GitHub page or email Tim Coorens (tcoorens@broadinstitute.org) or Mike Spencer Chapman (ms56@sanger.ac.uk).

### Installation

The necessary R package dependencies will be installed upon running the script. It is required to provide a path to an installation of [MPBoot](http://www.iqtree.org/mpboot/) to construct the phylogenetic tree topology. The path to MPBoot can be specified using the `--mpboot_path` parameter. This script has been verified to run on R3.6.1 and higher.

### Running the script
At the bare minimum, the script needs a matrix of total depth (variant sites by samples) and a matrix of number of variant-supporting reads (variant sites by samples). These can either be provided separately (using the `-r` and `-v` input parameters) or as a single output .tsv file from [cgpVAF](https://github.com/cancerit/vafCorrect), using the `-c` parameter. 

Using the example input files provided here, the command to run the script would be
```
Rscript build_phylogeny.R -r PD45567_NR.tsv -v PD45567_NV.tsv
```
Or
```
Rscript build_phylogeny.R -c PD45567.snp.tsv.gz
```

### Input options

Besides the necessary input files, all parameters and options can be adjusted as necessary. The full list of options can be accessed using `Rscript build_phylogeny.R -h`

```
Usage: build_phylogeny.R [options]


Options:
	-i DONOR_ID, --donor_id=DONOR_ID
		Patient/donor ID to add to names of output files

	-v INPUT_NV, --input_nv=INPUT_NV
		Input NV matrix (rows are variants, columns are samples)

	-r INPUT_NR, --input_nr=INPUT_NR
		Input NR matrix (rows are variants, columns are samples)

	-c CGPVAF_OUTPUT, --cgpvaf_output=CGPVAF_OUTPUT
		CGPVaf output file, instead of NR/NV matrices - can be multiple files, i.e. indel and snv data for the same donor (comma-separated)

	-o OUTPUT_DIR, --output_dir=OUTPUT_DIR
		Output directory for files

	-b BETA_BINOM_SHARED, --beta_binom_shared=BETA_BINOM_SHARED
		Only run beta-binomial filter on shared mutations. If FALSE, run on all mutations, before germline/depth filtering

	-n NCORES, --ncores=NCORES
		Number of cores to use for the beta-binomial step

	--normal_flt=NORMAL_FLT
		Name of the dummy normal to exclude from cgpVAF output

	--snv_rho=SNV_RHO
		Rho value threshold for SNVs

	--indel_rho=INDEL_RHO
		Rho value threshold for indels

	--min_cov=MIN_COV
		Lower threshold for mean coverage across variant site

	--max_cov=MAX_COV
		Upper threshold for mean coverage across variant site

	--only_snvs=ONLY_SNVS
		If indel file is provided, only use SNVs to construct the tree (indels will still be mapped to branches)

	--split_trees=SPLIT_TREES
		If both indels and SNVs are provided, plot trees separately for each.

	--keep_ancestral=KEEP_ANCESTRAL
		Keep an ancestral branch in the phylogeny for mutation mapping

	-x EXCLUDE_SAMPLES, --exclude_samples=EXCLUDE_SAMPLES
		Option to manually exclude certain samples from the analysis, separate with a comma

	--cnv_samples=CNV_SAMPLES
		Samples with CNVs, exclude from germline/depth-based filtering, separate with a comma

	--vaf_absent=VAF_ABSENT
		VAF threshold (autosomal) below which a variant is absent

	--vaf_present=VAF_PRESENT
		VAF threshold (autosomal) above which a variant is present

	-m MIXMODEL, --mixmodel=MIXMODEL
		Use a binomial mixture model to filter out non-clonal samples?

	-t TREE_MUT_PVAL, --tree_mut_pval=TREE_MUT_PVAL
		Pval threshold for treemut's mutation assignment

	-g GENOTYPE_CONV_PROB, --genotype_conv_prob=GENOTYPE_CONV_PROB
		Use a binomial mixture model to filter out non-clonal samples?

	-p MIN_PVAL_FOR_TRUE_SOMATIC, --min_pval_for_true_somatic=MIN_PVAL_FOR_TRUE_SOMATIC
		Pval threshold for somatic presence if generating a probabilistic genotype matrix

	--min_variant_reads_shared=MIN_VARIANT_READS_SHARED
		Minimum variant reads used in generating a probabilistic genotype matrix

	--min_vaf_shared=MIN_VAF_SHARED
		Minimum VAF used in generating a probabilistic genotype matrix

	--create_multi_tree=CREATE_MULTI_TREE
		Convert dichotomous tree from MPBoot to polytomous tree

	--mpboot_path=MPBOOT_PATH
		Path to MPBoot executable

	--germline_cutoff=GERMLINE_CUTOFF
		Log10 of germline qval cutoff

	--genomeFile=GENOMEFILE
		Reference genome fasta for plotting mutational spectra

	--plot_spectra=PLOT_SPECTRA
		Plot mutational spectra?

	--max_muts_plot=MAX_MUTS_PLOT
		Maximum number of SNVs to plot in mutational spectra

	--lowVAF_filter=LOWVAF_FILTER
		Minimum VAF threshold to filter out subclonal variants. Disabled by default.

	--lowVAF_filter_positive_samples=LOWVAF_FILTER_POSITIVE_SAMPLES
		Read number to apply exact binomial filter for samples with more than given number of reads. Disabled by default.

	--VAF_treshold_mixmodel=VAF_TRESHOLD_MIXMODEL
		VAF threshold for the mixture modelling step to consider a sample clonal

	-h, --help
		Show this help message and exit
```
### Examples

This command reproduces the results presented in the accompanying paper. The output files can also be found in the 'output' folder. Note that in order to use the "plot_spectra" feature, a path to the reference genome needs to be provided (`--genomeFile`). For the listed input, this is hg19/GRCh37. 

```
Rscript build_phylogeny.R -c PD45567.snp.tsv.gz -i PD45567 --exclude_samples PD45637b,PD45567f -m T --plot_spectra T --max_muts_plot 100000
```
