# Welcome to Metabolomic Selection for Enhanced Fruit Flavor

![icons of blueberries and a tomato](./fruit_icons.svg)

This repository contains data and scripts used to repoduce analyses in the manuscript "Metabolomic Selection for Enhanced Fruit Flavor" found on [BioRxiv](https://www.biorxiv.org/content/10.1101/2020.09.17.302802v1.full)

# Table of Contents
1. [Abstract](#Abstract)
2. [Figures](#figures)
3. [Conclusion](#conclusion)

# Introduction <a name="introduction"></a>

# Figures <a name="figures"></a>

Here we will go through the figures and which scripts were used to generate the underlying analysis. Often we generate the analysis in one script and design the figure component in another. We then combine the figure components together in inkscape.


## Figure 1 <a name="fig1"></a>

To generate this figure we first preprocess the input data, then create a network using WGCNA, then

* [0.preprocessing.R]   
We start by preprocessing the data from the supplemental files with default choices for imputation and scaling

* [1.a.wgcna_tomato.R]  
Next we create the metabolite network using the WGCNA package

* [1.b.metabolite_histograms.R]
Plotting the violin plots in panel b

![fig1](./figures/svgs/fig1.svg)

## Figure 2 <a name="fig2"></a>

![fig2](./figures/svgs/fig2.svg)

## Figure 3 <a name="fig3"></a>

![fig3](./figures/svgs/fig3.svg)

## Figure 4 <a name="fig4"></a>

![fig4](./figures/svgs/fig4.svg)

## Figure 5 <a name="fig5"></a>

![fig5](./figures/svgs/fig5.svg)

## Figure 6 <a name="fig6"></a>

![fig6](./figures/svgs/fig6.svg)





2.a.wgcna_blueberry.R
2.b.metabolite_histograms.R
3.a.variance_decomposition.R
4.a.1.metabolomic_selection_tomato.R
4.b.1.genomic_selection_tomato.sh
4.b.2.genomic_selection_tomato.sh
4.b.3.genomic_selection_tomato.R
4.b.4.metabolomic_selection_tomato.sh
4.b.5.metabolomic_selection_tomato.sh
4.b.6.metabolomic_selection_tomato.R
4.b.7.gblup_plots.R
4.c.1.subsampling.sh
4.c.2.subsampling.sh
4.c.3.subsampling.R
4.c.4.subsamplingPlots.R
5.a.1.calculate_final_weights.Rmd
5.a.2.plot_tomato_weights.R
5.a.3.plot_blueberry_weights.R
