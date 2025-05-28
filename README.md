# MRD-Adaptis
## ABSTRACT
<div style="text-align: justify;">
Detection of structural variations (SV) through circulating tumour DNA (ctDNA) has become a key method for detecting minimal residual disease (MRD). But The heterogeneity of ctDNA samples, characterized by limits of detection (LOD) and diverse structural variant (SV) types, causes inter- and intra-sample detection performance fluctuations, posing significant challenges for stable MRD detecting. SV detection tools are equipped with multiple user-configurable parameters designed to accommodate this heterogeneity, but optimizing these parameters is difficult for users due to the combinatorial explosion of possible parameter combinations and the sheer volume heterogeneous sequencing data. To address this, we propose MRD-Adaptis, a novel SV detection tool that incorporates a self-adaptive mechanism, optimizing parameters based on features derived from ctDNA sequencing data. This method integrates a Bayesian optimization framework with a meta-learning approach to dynamically adjust parameters, effectively adapting to sample heterogeneity and overcoming the complexities of manual parameter tuning. Extensive validation experiments using both simulated and real-world ctDNA datasets demonstrated that MRD-Adaptis consistently delivered high detection performance, as evidenced by substantially improved average F1-scores compared to existing methods. Additionally, MRD-Adaptis exhibited superior stability across heterogeneous ctDNA regions, demonstrated by reduced variance, lower root means square error (RMSE), and increased kurtosis. These results highlight the significant advantages of our self-adaptive parameter optimization framework in addressing sample heterogeneity, underscoring its potential to improve the accuracy and consistency of MRD detecting through ctDNA analysis.
</div>

<img src="https://github.com/aAT0047/MRD-Adaptis/blob/main/image/figure1.png" alt="Figure 1" width="600">


# Usage Guide

This guide explains the process of simulating data, training a meta-model, and testing models with recommended parameters using Python scripts and tools.

## Simulated Data Files Generation

Simulate fq data using **GSDcreator** by following the method described in the paper available on IEEE Xplore.

### Reference
[Simulation method for fq data](https://ieeexplore.ieee.org/abstract/document/8983192)

## Requirements
- **Python version:** 2.7 (for simulation steps)
- **Python version:** 3.6 or higher (for further data processing and meta-model tasks)

## Generate Simulated Data

### Generate Simulation Scripts
Use `A_stableCallerPaperSimFlow.py` to generate `.sh` scripts for simulating 10,000 samples:

    python A_stableCallerPaperSimFlow.py

### Install GSDcreator and Run Simulations
1. Convert scripts to Unix format and make them executable:

    dos2unix 10000run.py
    chmod +x 10000run.py

2. Execute the simulation script (ensure Python 2.7 is installed):

    python 10000run.py

## MRD-Adaptis Model
### step 1 Split Samples into Segments

Split samples into segments ranging from thousands (kilobases) to millions (megabases) of base pairs using Python 3.6+ with multithreading support:

    python aAT0047/MRD-Adaptis/blob/main/MRD-Adaptis/splitcsvbam.py

<img src="https://github.com/aAT0047/MRD-Adaptis/raw/main/image/figure2.png" alt="Figure 2: Sample Segmentation Process" width="600">

### step 2 Bayesian optimization framework 
#### Initializing opt Parameter Configuration (Multi-Threading)

Run the following script for initializing parameter configuration (Python 3.6+):

    python aAT0047/MRD-Adaptis/main.py

<img src="https://github.com/aAT0047/MRD-Adaptis/raw/main/image/figure3.png" alt="Figure 3: Meta-Model Training Process" width="600">

#### Extracting Sample Meta-Features

Use Python 3.6+ to extract sample meta-features:

    python aAT0047/MRD-Adaptis/metafeature.py
    
### step 3 Training a Meta-Model

Training a meta-model involves creating a model that learns from the outputs or performance of other models. The resulting meta-model is saved as `multi_target_regression_model.pth`.
## Meta-features

| Meta-feature                                    | Description                                                                          |
|-------------------------------------------------|--------------------------------------------------------------------------------------|
| The average length of reads                     | The mean length of sequencing reads in a given sub-region.                           |
| The average sequencing depth                    | The average coverage of sequencing reads in a given sub-region, reflecting the overall read density. |
| The average invert size                         | The mean size of inversions detected within the sub-region.                          |
| Average Mapping Quality                         | The average quality score of read alignments in the sub-region.                      |
| The proportion of short breakpoints             | The proportion of structural variations with short breakpoints in the sub-region (< 500). |
| The proportion of middle breakpoints            | The proportion of structural variations with medium-sized breakpoints in the sub-region (500–2000). |
| The proportion of large breakpoints             | The proportion of structural variations with large breakpoints in the sub-region (> 2000). |
| The proportion of breakpoints in repeat regions | The proportion of breakpoints located within repetitive genomic regions in the sub-region. |
| Average SV length                               | The average length of structural variations detected in the sub-region.              |
| Insert Size Mean                                | The mean insert size (distance between paired-end reads) within the sub-region.      |
| Insert Size Standard Deviation                  | The standard deviation of insert sizes, reflecting variability in paired-end distances. |
| Read Mapping Bias (RMB)                         | The proportion of reads that align asymmetrically across the sub-region, indicating potential biases. |
| HDMP                                            | High Mean Density Percentage (Number of Small Gaps / Total Number of Gaps).          |

## multi_target :Parameters

| Parameter Name              | Meaning                                                                                                                                         |
|-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `-w`                        | Window size: Defines the range of the detection window. Larger windows cover a broader range but may reduce precision.                          |
| `-msw`                      | Minimum Per-Sample Weight: The minimum weight required for each sample to call a variant.                                                        |
| `-tt`                       | Trim Threshold: Filters out overly long or short reads to ensure only qualified reads are processed, improving detection accuracy.              |
| `-back_distance`            | Backward Distance: The maximum distance the algorithm can check backward.                                                                        |
| `-min_mapping_threshold`    | Minimum Mapping Quality Threshold: Filters out low-quality mappings.                                                                             |
| `-min_clip`                 | Minimum Clipping Length: The minimum length of clipped segments to be considered potentially related to structural variants.                     |
| `-read_length`              | Read Length: The length of each read in genome sequencing. Used for normalizing sequencing data.                                                 |
| `-min_non_overlap`          | Minimum Non-Overlap: The minimum non-overlapping region between two sequences required to identify structural variants.                         |
| `-discordant_z`             | Discordant Z-Score: Measures whether paired reads significantly deviate from expected positions or orientation.                                 |


    python aAT0047/MRD-Adaptis/metaleaner.py

<img src="https://github.com/aAT0047/MRD-Adaptis/raw/main/image/figure4.png" alt="Figure 4: Meta-Model Evaluation" width="600">

## Testing Model with Recommended Parameters

Run the following command to test the model with recommended parameters:

    python model.py 1.bam

**Example Output:**  
The `prediction_dict` is generated with keys such as:

| Parameter                 | Value |
|---------------------------|-------|
| `w`                       | …     |
| `msw`                     | …     |
| `tt`                      | …     |
| `back_distance`           | …     |
| `min_mapping_threshold`   | …     |
| `min_clip`                | …     |
| `read_length`             | …     |
| `min_non_overlap`         | …     |
| `discordant_z`            | …     |


## Testing DELLY, LUMPY, Manta, BreakDancer, Pindel, MetaSV, SvABA & MRD-Adaptis

To evaluate the model's generalization ability, we trained it on 400 simulated data samples and tested it using 176 real-world data samples. Additionally, 5-fold cross-validation was performed on the training set to ensure robustness and prevent overfitting.Use the following script to test different tools:

    python /SVfolder/vsworkflow/callerworkflow.py
<img src="https://github.com/aAT0047/MRD-Adaptis/raw/main/image/figure7.png" alt="Figure 7: Prediction Example" width="600">

<img src="https://github.com/aAT0047/MRD-Adaptis/raw/main/image/All_Metrics_ismb.png"  width="600">
