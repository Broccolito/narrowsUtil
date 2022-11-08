# NarrowsUtil

NarrowsUtil - Utilities for the Narrows cluster of the Salem Lab



## Installation

In the narrows cluster, create a conda environment using:

```bash
# Update the conda base environment to stay up-to-date
conda update -n base -c defaults conda

# Check if the narrowsUtil environment has already been created
conda env list

# If not, create the narrowsUtil conda environment
conda create --name narrowsUtil

# Activate the environment
conda activate narrowsUtil

# Update the conda environment
conda update --all

# Install R packages and other dependencies needed for narrowsUtil
conda install -c conda-forge python=3.10 -y
conda install -c conda-forge r-base -y
conda install -c anaconda cmake -y
conda install -c conda-forge r-devtools -y
conda install -c conda-forge r-ggplot2 -y
conda install -c conda-forge r-ggpubr -y
conda install -c conda-forge r-data.table -y
```



## Usage

#### Merge SAIGE output

```bash
# Merge chromosomal summary stats
Rscript /salemlab/users/wagu/narrowsUtil/merge_saige.R \
	HAPO_AS_increment_chrXXX.txt TRUE
	
# With chrX results
Rscript /salemlab/users/wagu/narrowsUtil/merge_saige.R \
	HAPO_AS_increment_chrXXX.txt TRUE \
	HAPO_AS_increment_plink2_run_chrX.txt TRUE
```

#### Clumping

```bash
# Clumping the GWAS signals to SAS
bash /salemlab/REF/SOFTWARE/SHCode/Filter_and_Clump_GWAS_v1.sh \
  HAPO_AS_increment_chr_MERGED.txt \
  P_gc SNPID SAS \
  HAPO_AS_increment_CLUMPED
```

#### Annotate GWAS summary statistics

```bash
# Example
Rscript /salemlab/users/wagu/narrowsUtil/annotate_gwas.R \
  HAPO_MA_increment_CLUMPED_1e5_snplist.txt \
  HAPO_MA_increment_annotated.xlsx
```

#### Make Manhattan plot and Q-Q plot

```bash
# Example
Rscript /salemlab/users/wagu/narrowsUtil/make_manhqq_plot.R \
  HAPO_AS_increment_chr_MERGED.txt \
  HAPO_AS_increment
```

#### Make Forestplot

```bash
# Arguments
Rscript /salemlab/users/wagu/narrowsUtil/make_senlin_plot.R \
	[Args1 METAL Parameter File] \
	[Args2 Meta Analysis Result] \
	[Args3 SNPID] \
	[Args4 plot and stats name suffix]

# Example
Rscript /salemlab/users/wagu/narrowsUtil/make_senlin_plot.R \
  /salemlab/users/m1ma/forestPlot/testsenlin.par \
  /salemlab/users/wagu/hapo/minimal_model/HAPO_BINARY_META1.txt \
  6:150541053:G:T \
  forestplot_test
```

#### Make Forestplot from customized stats sheet

```R
Rscript /salemlab/users/wagu/narrowsUtil/make_senlin_plot_from_stats.R \
	[Args1 senlin plot parameter file] \
	[Args2 plot name]

# Example
Rscript /salemlab/users/wagu/narrowsUtil/make_senlin_plot_from_stats.R \
	/salemlab/users/wagu/narrowsUtil/forestplot_test_forestplot_stats.csv \
	example_plot.png
```

#### Make Crosstrait Forestplot

```r
# Arguments
Rscript /salemlab/users/wagu/narrowsUtil/make_crosstrait_senlin_plot.R \
	[Args1 METAL Parameter File] \
	[Args3 SNPID] \
	[Args4 plot and stats name suffix]

# Example
Rscript /salemlab/users/wagu/narrowsUtil/make_crosstrait_senlin_plot.R \
  /salemlab/users/m1ma/forestPlot/testcross.par \
  6:150541053:G:T \
  crosstraitplot_test
```



#### Make Crosstrait Forestplot from customized stats sheet

```r
Rscript /salemlab/users/wagu/narrowsUtil/make_crosstrait_senlin_plot_from_stats.R \
	[Args1 senlin plot parameter file] \
	[Args2 plot name]

# Example
# Need to check output file name from last part to be arguments for following step
Rscript /salemlab/users/wagu/narrowsUtil/make_crosstrait_senlin_plot_from_stats.R \
	/salemlab/users/wagu/narrowsUtil/crosstraitplot_test_forestplot_stats.csv \
	example_plot.png
```





### Specifying Cluster Parameter

```bash
#! /bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name  HAPO_INCREMENT_RUNNER
#SBATCH --output=HAPO_INCREMENT_RUNNER.out
#SBATCH --partition=salem-compute

# Activate conda environment
conda activate narrowsUtil

Rscript /salemlab/users/wagu/narrowsUtil/make_senlin_plot.R testsenlin.par /salemlab/users/wagu/toy_meta/HAPO_BINARY_META1.txt 6:150541053:G:T forestplot_test



