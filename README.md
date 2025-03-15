# pi-USEA: Ubiquitin Ligase/Deubiquitylase Substrate Enrichment Analysis

This tool is specifically designed to predict the activity of ubiquitin ligases and deubiquitylases. Users can upload the ratio change data of the ubiquitination group of interest to obtain corresponding enzyme activity prediction results.

## Features

1.Dual Activity Prediction: Supports the prediction of ubiquitination (E3 ligase) and deubiquitination (DUBs) activities.
2.Flexible Input: Accepts data with or without site-specific modification information.
3.Dual Prediction Mode: Supports analysis based on validated interaction data, and can also integrate predicted interaction data for combined calculations.
4.Multi-Species Support: Pre-configured datasets for human and mouse systems.
5.Validated Performance: Rigorously tested on benchmark datasets, ensuring reliable results.

## Usage

1.Environment Setup: It is recommended that users download the code and dataset files to their local R environment. Code download link: [pi-USEA.R](https://github.com/PHOENIXcenter/pi-USEA/blob/main/code/pi-USEA.R), Dataset download link: [Dataset](https://github.com/PHOENIXcenter/pi-USEA/tree/main/dataset).

<div style="text-align: center;">

# Dataset Sources

![alt text](image-1.png)

</div>

2. Data Format Requirements: Please modify your ubiquitination group ratio change data table to the `test_data.csv` format, as shown below:

   | SUB_ACC_ID | SYMBOL | SUB_MOD_RSD | p          | FC            |
   |------------|--------|--------------|------------|---------------|
   | A2A432     | Cul4b  | K402         | 0.411203714 | -4.376820087  |
   | A2A432     | Cul4b  | K460         | 0.179046853 | -10.82791042  |
   | A2A432     | Cul4b  | K495         | 0.558549777 | -0.19614013   |
   | A2A6A1     | Gpatch8| K1341       | 0.237466694 | 3.840063575   |
   | A2A6Q5     | Cdc27  | K133         | 0.000778651 | -10.60901848  |

   Example file: [test_data.csv Example](https://github.com/PHOENIXcenter/pi-USEA/blob/main/test/test_data.csv)

## Notes

- This tool's dataset is organized specifically for human and mouse species. If your data is from human or mouse, you can directly select the corresponding species to download and use the built-in dataset for analysis. If it is mouse data, you can skip the data integration step and directly download the `part3_calculation` dataset ([Download Link](https://github.com/PHOENIXcenter/pi-USEA/tree/main/dataset/part3_calculation)), and start running from `part3` after installing the required R packages.
  
- If you want to analyze species other than human and mouse, please download the relevant data yourself, and refer to the format of this tool's dataset to modify the column names and save it as a `.csv` file. Recommended database download links:
  - [UbiBrowser](http://ubibrowser.bio-it.cn/ubibrowser_v3/home/download)
  - [UbiNet](https://awi.cuhk.edu.cn/~ubinet/download.php)
  - [E3Net](http://pnet.kaist.ac.kr/e3net)

## Interaction Dataset Format Recommendations

Please ensure that the predicted interaction dataset contains the enzyme protein ID (`ENZY_ACC_ID`) and the substrate protein ID (`SUB_ACC_ID`), formatted as follows:

| ENZY_ACC_ID | SUB_ACC_ID |
|-------------|------------|
| P46737      | Q3UPL0     |
| P46737      | Q6ZQA6     |
| P46737      | E0CZ16     |
| P46737      | Q80XI4     |
| P46737      | Q8K2P1     |

The real interaction dataset should include the enzyme protein ID (`ENZY_ACC_ID`), substrate protein ID (`SUB_ACC_ID`), and substrate gene ID (`SUB_GENE`), formatted as follows:

| ENZYME | ENZY_ACC_ID | GENE | KIN_ORGANISM | SUBSTRATE | SUB_GENE_ID | SUB_ACC_ID | SUB_GENE | SUB_ORGANISM |
|--------|-------------|------|--------------|------------|--------------|------------|----------|--------------|
| AHR    | P30561      | Ahr  | mouse        | PPARG      | P37238       | Pparg     | mouse        |
| AHR    | P30561      | Ahr  | mouse        | RIPK1      | Q60855       | Ripk1 | mouse        |
| AMFR   | Q9R049      | Amfr | mouse        | CD82       | P40237       | Cd82      | mouse        |
| ARI1   | Q9Z1K5      | Arih1| mouse        | PD1L1      | Q9EP73       | Cd274     | mouse        |
| ASB17  | Q8VHP9      | Asb17| mouse        | MCL1       | P97287       | Mcl1      | mouse        |

## Site Dataset Format Recommendations

Please ensure that the site dataset contains all ubiquitination sites, formatted as follows:

### PSP_site.csv
| ENZYME | ENZY_ACC_ID | SUB_MOD_RSD |
|--------|-------------|--------------|
| Ywhab  | Q9CQV8      | K11          |
| Ywhab  | Q9CQV8      | K13          |
| Ywhab  | Q9CQV8      | K29          |
| Ywhab  | Q9CQV8      | K51          |

### Ubisite.csv
| Uniprot ID | Sites |
|------------|-------|
| O00255     | 120   |
| O00255     | 315   |
| O00255     | 4     |
| O00255     | 609   |

### IUUCD.csv
| Gene_name | Ubiquitination_site |
|-----------|---------------------|
| Rag1      | 233                 |
| Zranb3    | 322                 |
| Psmd4     | 135;122;74;40;126;364 |

### TransDSI.csv
| SwissProt ID (DUB) | SwissProt ID (SUB) | Entrez Gene Symbol (DUB) | Entrez Gene Symbol (SUB) | TransDSI Score | Contributing Residues (DUB) | Contributing Residues (SUB) |
|---------------------|---------------------|---------------------------|---------------------------|----------------|------------------------------|------------------------------|
| A6NNY8              | Q9HAW4              | USP27X                    | CLSPN                     | 0.992          | [S71, S72, F73, T74, P116, G159, T262, T325, Y326, I327] | [K179, E181, K743, K930, D932, R1040, E1042, K1321, T1322, D1323] |
| A6NNY8              | P20226              | USP27X                    | TBP                       | 0.969          | [I101, S128, L129, F130, T336, S365, L366, F367, S383, S384] | [Q67, Q74, Q75, Q76, Q77, Q78, Q79, Q80, Q82, Q87] |
| A6NNY8              | P62380              | USP27X                    | TBPL1                     | 0.969          | [F73, T74, I75, P116, S117, D200, T262, P303, I327, L423] | [I11, L12, I13, K30, V44, R82, R83, I96, K101, H135] |
| A6NNY8              | Q13547              | USP27X                    | HDAC1                     | 0.968          | [S71, S72, F73, P116, T221, T262, L263, T325, Y326, I327] | [K50, E52, M249, M251, S263, S421, D422, S423, T445, E446] |

## Contact Information

If you have any questions during usage, please feel free to contact us:

- Chang Cheng: changchengbio@163.com or changcheng@ncpsb.org.cn
- Liu Ning: liu_ning2021@163.com
