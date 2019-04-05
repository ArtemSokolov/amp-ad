# Informatics analysis of AMP-AD datasets to identify drug repurposing candidates for Alzheimers Disease

Laboratory of Systems Pharmacology, Harvard Medical School

The computational pipeline contained in this repository works with several sources of data. The first is mRNA expression collected *in vitro* from differentiated human neural progenitor cells treated with a panel of compounds. Some of these compounds are FDA-approved for non-AD indications, other compounds are in various stages of a clinical trial. The data is processed by a standard differential gene expression tool ([edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)) to arrive at a set of genes associated with each drug.

The resulting gene sets are profiled against the second source of data, [AMP-AD](http://www.synapse.org/AMPAD), which contains publicly available -omics data collected on postmortem brain tissue of Alzheimer's patients. The AMP-AD data is used to build predictors of disease stage from the expression of genes identified in the *in vitro* experiments. The predictors are systematically evaluated across several AMP-AD datasets, brain regions, and binary classification tasks using the [BTR](https://github.com/pvtodorov/btr) framework on a local compute cluster. The resulting analysis produces a list of drugs ranked by the strength of association between their molecular mechanisms of action and disease. Strong associations may be due to potential neuroprotective or neurotoxic effects, and efforts to experimentally identify the direction of toxicity are underway.

The third part of the analysis utilizes target affinity spectrum data from [The Small Molecule Suite](https://labsyspharm.shinyapps.io/smallmoleculesuite/) to ask what combinations of targets drive the position of drugs in the ranked list derived above. By considering both FDA-approved and non-approved drugs, we expand the space of molecular mechanisms, allowing us to generate novel hypotheses about actionable biology of the disease. These hypotheses will guide our selection of the compound panel to profile next. In parallel, we are working on designing experiments to confirm/reject the neuroprotective effect of top-ranking compounds. The latest results are described in our manuscript, which is currently in preparation.

## Funding

We gratefully acknowledge support by NIA grant 1 R01 AG058063-01A1: Harnessing Diverse BioInformatic Approaches to Repurpose Drugs for Alzheimers Disease
