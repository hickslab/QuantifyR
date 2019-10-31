# `QuantifyR`

## Instructions

1. [Read the `QuantifyR` tutorial](https://hickslab.github.io/QuantifyR/)
2. Collect the required data.
3. Download a [workflow](https://github.com/hickslab/QuantifyR/tree/master/workflow) based on the experiment.
4. Find and update sections marked by `???` in the workflow.

## Required

### Progenesis:

|Experiment|Peptide Measurements|Protein Measurements|Database|
|:-:|:-:|:-:|:-:|
|Global||[X](https://raw.githubusercontent.com/hickslab/QuantifyR/master/data/20180502_WOS52_Cr_UPS_protm.csv)||
|PTM|[X](https://raw.githubusercontent.com/hickslab/QuantifyR/master/data/20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_pepm.csv)|[X](https://raw.githubusercontent.com/hickslab/QuantifyR/master/data/20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_protm.csv)|[X](https://raw.githubusercontent.com/hickslab/QuantifyR/master/data/Cr_uniprot_crap_20190130.fasta)|

## Workflow

### LFQ:

#### Progenesis QI for proteomics v2.0
* [Global](https://raw.githubusercontent.com/hickslab/QuantifyR/master/workflow/Global-LFQ.R)
* [Peptide](https://raw.githubusercontent.com/hickslab/QuantifyR/master/workflow/Peptide-LFQ.R)
* [Ox-RAC](https://raw.githubusercontent.com/hickslab/QuantifyR/master/workflow/OxRAC-LFQ.R)
* [Phosphorylation](https://raw.githubusercontent.com/hickslab/QuantifyR/master/workflow/Phospho-LFQ.R)

