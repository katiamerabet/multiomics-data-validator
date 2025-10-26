(meta_data_env) katiamerabet@katiamerabet-latitudee7450:~/Downloads/Meta_data_validation_tool$ cat README.md
# Multi-Omics Data Validator

An interactive **Streamlit web app** that validates the **completeness, consistency, and structure** of multi-omics datasets (genomics, proteomics, metabolomics, etc.).

This tool automatically checks your data for:
- Missing values  
- Duplicates  
- Outliers  
- Data consistency issues (e.g., negative values)  
- Categorical variable validation  
- Distributions and visualizations  

---

##  Motivation

This project was developed as part of my ongoing interest in the intersection of AI, data science, and bioinformatics, with the goal of improving the quality and reproducibility of biological research data. By building a user-friendly validation tool with Streamlit, I wanted to facilitate better metadata integrity, enabling researchers to trust and analyze their datasets more efficiently. This project also reflects my broader motivation to apply machine learning and data-driven methods to address real challenges in life sciences and research environments.

---

##  Features

- Upload your dataset in `.csv`, `.tsv`, or `.xlsx` format.  
- Automatically detects issues (missing data, duplicates, outliers).  
- Generates a **validation report** summarizing dataset quality.  
- Displays **visualizations** of numeric columns.  
- Ideal for **multi-omics quality control** and **data preprocessing** before analysis.

---

##  Example Use Case

This app can be used to validate data before:
- Multi-omics integration (e.g., transcriptomics + proteomics)
- Machine learning pipelines
- Statistical or differential expression analyses
- Data submission or publication QC

---

# Screenshots

Here are some screenshots:

<p align="center">
  <img src="assets/Screenshot_20251026_111254.png" width="45%"/>
  <img src="assets/Screenshot_20251026_111543.png" width="45%"/>
</p>
---

##  Installation

Clone the repository:
```bash
git clone https://github.com/katiamerabet/multiomics-data-validator.git
cd multiomics-data-validator