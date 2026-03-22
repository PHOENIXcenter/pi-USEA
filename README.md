# pi-USEA: Ubiquitin Ligase/Deubiquitase Substrate Enrichment Analysis

## 1. Tool Overview

This tool is specifically designed to predict the activity of ubiquitin ligases and deubiquitinases. Users only need to upload the ratio change data of their ubiquitinome of interest to obtain corresponding enzyme activity prediction results.

## 2. Core Features

1. **Dual Activity Prediction**: Supports prediction of both ubiquitin (E3 ligase) and deubiquitinase (DUB) activity.
2. **Flexible Input**: Accepts data with or without site-specific modification information.
3. **Dual Prediction Modes**: Supports analysis based on validated interaction data, and can also integrate predicted interaction data for weighted combined calculations.
4. **Multi-Species Support**: Pre-loaded reference datasets for human (*Homo sapiens*) and mouse (*Mus musculus*).
5. **Validated Performance**: Rigorously tested on benchmark datasets, ensuring reliable results.

## 3. Usage Workflow

### 3.1 Environment Setup

- **R Environment**: Requires R version 3.5 or higher.
- **Code and Dataset Download**:
    - Code files: [[Address](https://github.com/PHOENIXcenter/pi-USEA/blob/main/code/pi-USEA)]
    - Datasets: part1_interaction[[Address](https://github.com/PHOENIXcenter/pi-USEA/tree/main/dataset/part1_interaction)], part2_calculation:[ Baidu Netdisk](https://pan.baidu.com/s/1KQOuvBGNUoS0UCzCUuk4TA?pwd=mp94) (Extraction Code: mp94)

### 3.2 Input Data Preparation

- **File Format**: Can be saved as `test_data.csv`, containing the following columns (at minimum, `SUB_ACC_ID` and `FC` are required):

| Column Name | Description | Example |
| :--- | :--- | :--- |
| SUB_ACC_ID | UniProt ID of the substrate protein | A2A432 |
| SYMBOL | Substrate gene symbol (optional) | Cul4b |
| SUB_MOD_RSD | Substrate modification site (optional) | K402p |
| P.Value | Adjusted p-value (optional) | 0.411203714 |
| FC | log2 ratio of experimental group/control group (required) | -4.376820087 |

   - Example file: [test_data.csv Example](https://github.com/PHOENIXcenter/pi-USEA/blob/main/test/test_data.csv)

### 3.3 Analysis Steps

The tool's datasets have been pre-processed for human and mouse species. If your data is from human or mouse, you can directly use the corresponding species dataset provided by the tool. Open `pi-USEA.R`, install and load the required R packages, and start the analysis directly from Part 2. The main steps are:

1.  **Read Interaction Datasets**: Load the interaction datasets (including `E3pre`: combined predicted and real E3 interactions; `E3real`: real E3 interactions; `DUBpre`: combined predicted and real DUB interactions; `DUBreal`: real DUB interactions) and your prepared ubiquitinome data.
2.  **Merge Datasets**: Merge the interaction datasets with your ubiquitinome data.
3.  **Analyze Real Datasets (`E3real`, `DUBreal`)**: Directly calculate enzyme activity and output results `E3/DUB_real_activity.csv`.
4.  **Analyze Predicted Datasets (`E3pre`, `DUBpre`)**: Perform confidence assessment, select a confidence threshold, and filter high-confidence interaction entries (manual threshold setting required) for subsequent analysis, for example:
    `E3pre <- E3pre %>% filter(interScore > 0.75)`
5.  **Weighted Calculation for Predicted Datasets**: Calculate enzyme activity using a weighted average and output results.
    - `E3_real_activity.csv` / `DUB_real_activity.csv`: Activity results based on validated interactions.
    - `E3_pre_activity.csv` / `DUB_pre_activity.csv`: Activity results integrating predicted data.

## 4. Result Field Descriptions

### 4.1 Core Identifiers

| Field Name | Description | Example |
| :--- | :--- | :--- |
| ENZY_ACC_ID | UniProt ID of the enzyme (E3 ubiquitin ligase or DUB) | P04637 (human p53 protein) |
| SUB_ACC_ID | UniProt ID of the substrate protein | Q13547 (histone deacetylase 1) |
| SUB_MOD_RSD.x | Modified residue position on the substrate protein (merged from `SUB_MOD_RSD` in test_data) | K48, K341, K363 |
| ENZYME | Name of the enzyme | UBE3A (E3 ubiquitin ligase) |
| SUBSTRATE | Name of the substrate protein | TP53 (tumor suppressor protein p53) |
| GENE | Gene name corresponding to the enzyme | Ube3a |
| SUB_GENE | Gene name corresponding to the substrate | Hdac1 |
| SYMBOL | Gene symbol (from test_data) | Hdac1 |

### 4.2 Experimental Data

| Field Name | Description | Calculation |
| :--- | :--- | :--- |
| FC | Fold change (from test_data) | log2 (experimental group / control group) |
| logFC_mean | Average log2FC for all substrate sites corresponding to a single enzyme | mean(log2(FC)) |
| logFC_count | Number of substrate sites corresponding to a single enzyme | Direct count |

### 4.3 Statistical Significance

| Field Name | Description | Formula |
| :--- | :--- | :--- |
| z_score | Quantitative result of enzyme activity. Positive indicates increased activity, negative indicates decreased activity. | Real: (logFC_mean - overall mean) / (overall sd / sqrt(n)) <br> Predicted: (weighted logFC_mean - overall weighted mean) / (overall weighted sd / sqrt(weight)) <br> Note: DUB direction is opposite. |
| p_value | p-value calculated from z_score | `2 * pnorm(-abs(z_score))` (two-tailed test) |
| fdr | False discovery rate corrected p-value | Benjamini-Hochberg correction |

### 4.4 Prediction Score Fields (This analysis primarily uses `interScore`. For details, refer to the UbiBrowser original publication.)

| Field Name | Description | Score Range | Source |
| :--- | :--- | :--- | :--- |
| interScore | Enzyme-substrate interaction prediction score | 0.5 - 1 | UbiBrowser and other prediction tools |
| domainScore | Domain interaction score | - | InterPro and other domain databases |
| motifScore | Substrate motif matching score | - | e.g., linear ubiquitination motif matching |
| netScore | Network topology importance score | - | Centrality analysis based on PPI networks |
| goScore | Gene ontology functional similarity score | - | GO semantic similarity calculation |

### 4.5 Source Markers

| Field Name | Description | Typical Values |
| :--- | :--- | :--- |
| source | Data source (after merging) | UbiBrowser_real, E3Net, UbiBrowser_pre |
| source_pre | Predicted data source | UbiBrowser_pre |
| source.x / source.y | Redundancy columns for real vs. predicted interactions | x: UbiBrowser_pre <br> y: UbiNet, UbiBrowser_real |

### 4.6 Threshold Evaluation

| Field Name | Description | Purpose |
| :--- | :--- | :--- |
| Threshold | `interScore` score threshold | Used to filter reliable interactions |
| num_unique_enzy | Number of unique enzymes at the current threshold | Assess coverage |
| num_rows | Total number of interactions at the current threshold | Assess data volume |

### 4.7 Other Technical Fields

| Field Name | Description |
| :--- | :--- |
| KIN_ORGANISM | Species of the enzyme |
| SUB_ORGANISM | Species of the substrate |
| exists_real | Indicates if the interaction is real (TRUE/FALSE) |
| normalized_interScore | Normalized interaction score (`(interScore - min) / (max - min)`) |

## 5. Custom Analysis (Non-Human/Mouse Species)

1.  If you need to analyze species other than human and mouse, please download interaction data from suggested databases:
    - [UbiBrowser](http://ubibrowser.bio-it.cn/ubibrowser_v3/home/download)
    - [UbiNet](http://csb.cse.yzu.edu.tw/UbiNet/)
    - [E3Net](http://pnet.kaist.ac.kr/e3net)

2.  Modify column names according to the format used in this tool's datasets:
    - **Real interaction datasets**: After downloading, modify column names to include `ENZY_ACC_ID`: UniProt ID of the enzyme (required), `SUB_ACC_ID`: UniProt ID of the substrate (required).
    - **Predicted interaction datasets from UbiBrowser**: Can be read directly from the downloaded `.gz` file.
    - After preparing the required file formats and names, install necessary packages and start the analysis from Part 1.


## Contact Information

If you have any questions during usage, please feel free to contact us:

- **Chang Cheng** : changchengbio@163.com or changcheng@ncpsb.org.cn  
- **Liu Ning**: liu_ning2021@163.com