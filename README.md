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



#### Annotate GWAS summary statistics



#### Make Manhattan plot and Q-Q plot



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
  /salemlab/users/wagu/hapo/minimal_model/meta_analysis.par \
  /salemlab/users/wagu/hapo/minimal_model/HAPO_BINARY_META1.txt \
  6:150541053:G:T \
  forestplot_test
```







