# Bioinformatic analysis of metabarcoding data for soil invertebrates

## 1. Prerequisites

### 1.1. Install required softwares

- UNIX environment
- R and Rstudio
- QIIME2
- FastQC (optional)
- MultiQC (optional)

#### 1.1.1. Obtain a Unix environment
A **Unix** environment is an operating system interface that allows users to manage files and run programs through a command line, it is widely used for scientific computing and data analysis. **Linux** and **MAC** operating systems are already Unix-based. For **Windows** users the easiest way to obtain a working Unix environment is to install the Windows Subsystem for Linux (WSL) following these [instructions](https://learn.microsoft.com/en-us/windows/wsl/install). Then to activate the Linux terminal you just need to type `wsl` in the windows PowerShell.

#### 1.1.2. Install R and RStudio
**R** is a programming language designed for statistical computing and data analysis, particularly well suited for handling, analyzing, and visualizing biological data. To install **R** follow the instructions for your operating system reported in this [page](https://cran.rstudio.com/).

**RStudio** is an integrated development environment (IDE) that provides a user-friendly interface for writing, running, and managing R code, along with tools for data visualization and debugging. It is not strictly needed to run the analysis but it makes it easier. You can install **RStudio** following the instructions for your operating system reported in this [page](https://posit.co/download/rstudio-desktop/).

#### 1.1.3. Install QIIME2
**QIIME2** is an open-source, plugin-based bioinformatics platform for reproducible analysis of amplicon sequencing data, widely used to study microbial and eukaryotic community composition from DNA metabarcoding datasets. The easiest way to install QIIME2 is within a conda environment, as described below, if you prefer to install it within a docker follow these [instructions](https://library.qiime2.org/quickstart/amplicon#:~:text=Using%20Docker%C2%B6,%2Dt) instead.

**STEP 1 - Install Miniconda**

Follow the instructions for your operating system. In case of trubles, the full istructions for miniconda installation are available [here](https://www.anaconda.com/docs/getting-started/miniconda/install#manual-shell-initialization).

- **Linux or Windows Subsystem for Linux**

    1. Download the latest Miniconda installer.

        ```bash
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        ```

    2. Verify installer integrity (this code is valid for the installer release 2025-12-16, check the sha256 value for different versions [here](https://repo.anaconda.com/miniconda/)). If the output of the command is 'OK' go to the next instruction. If the output is 'MISMATCH' something went wrong so you should delete the installer file with `rm Miniconda3-latest-Linux-x86_64.sh` and download it again with the previous command.
    
        ```bash
        [[ "$(sha256sum Miniconda3-latest-Linux-x86_64.sh | cut -d' ' -f1)" == "e0b10e050e8928e2eb9aad2c522ee3b5d31d30048b8a9997663a8a460d538cef" ]] && echo "OK" || echo "MISMATCH"
        ```

    3. Install Miniconda. During installation follow the instructions. Enter `yes` to accept licence agreement. Press Return to accept the default install location (`PREFIX=/Users/<USER>/miniconda3`), or enter another file path to specify an alternate installation directory. The installation might take a few minutes to complete. At the end of the installation process, enter `yes` to initialize conda so that whenever you open a new shell it recognizes conda commands automatically.
    
        ```bash
        bash Miniconda3-latest-Linux-x86_64.sh
        ```
    
    4. Delete the installer.
    
        ```bash
        rm Miniconda3-latest-Linux-x86_64.sh
        ```
    
    5. ⚠️ You need to close and re-open your terminal window for the installation to fully take effect.
    
    6. Verify the installation. If everything is ok it should print the version of conda installed.
    
        ```bash
        conda --version
        ```


- **MAC OS**

    1. Download the latest Miniconda installer.
    
        ```bash
        curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
        ```
    
    2. Verify installer integrity (this code is valid for the installer release 2025-12-16, check the sha256 value for different versions [here](https://repo.anaconda.com/miniconda/)). If the output of the command is 'OK' go to the next instruction. If the output is 'MISMATCH' something went wrong so you should delete the installer file with `rm Miniconda3-latest-Linux-x86_64.sh` and download it again with the previous command.
    
        ```bash
        [[ "$(shasum -a 256 Miniconda3-latest-MacOSX-arm64.sh | cut -d' ' -f1)" == "9f84ad10ea513fb59bb714933bc8dc092bd25fdb03c236868f5d5af3c26a1fd4" ]] && echo "OK" || echo "MISMATCH"
        ```
    
    3. Install Miniconda. During installation follow the instructions. Enter `yes` to accept licence agreement. Press Return to accept the default install location (`PREFIX=/Users/<USER>/miniconda3`), or enter another file path to specify an alternate installation directory. The installation might take a few minutes to complete. At the end of the installation process, enter `yes` to initialize conda so that whenever you open a new shell it recognizes conda commands automatically.
    
        ```bash
        bash Miniconda3-latest-MacOSX-arm64.sh
        ```
    
    4. Delete the installer.
    
        ```bash
        rm Miniconda3-latest-MacOSX-arm64.sh
        ```
    
    5. ⚠️ You need to close and re-open your terminal window for the installation to fully take effect.
    
    6. Verify the installation. If everything is ok it should print the version of conda installed.
    
        ```bash
        conda --version
        ```



**STEP 2 - Install QIIME2**

Follow the instructions for your operating system. 

- **Linux or Windows Subsystem for Linux**

    1. Update conda to the latest version.
    
        ```bash
        conda update conda
        ```
    
    2. Create a new conda environment for QIIME 2 and install it. You can choose whatever name you’d like for the environment, here we call it 'qiime2-amplicon-2026.1' to indicate what QIIME2 release is installed.
    
        ```bash
        conda env create --name qiime2-amplicon-2026.1 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
        ```
    
    3. Activate the conda environment containing QIIME2. You will need to run this command every time you need to use QIIME2.
    
        ```bash
        conda activate qiime2-amplicon-2026.1
        ```
    
    4. Deactivate the conda environment. Run this to quit the QIIME2 environment and return to the base environment.
    
        ```bash
        conda deactivate
        ```

- **MAC OS (Intel)**

    1. Update conda to the latest version.
    
        ```bash
        conda update conda
        ```
    
    2. Create a new conda environment for QIIME 2 and install it. You can choose whatever name you’d like for the environment, here we call it 'qiime2-amplicon-2025.4' to indicate what QIIME2 release is installed.
    
        ```bash
        conda env create --name qiime2-amplicon-2026.1 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
        ```
    
    3. Activate the conda environment containing QIIME2. You will need to run this command every time you need to use QIIME2.
    
        ```bash
        conda activate qiime2-amplicon-2026.1
        ```
    
    4. Deactivate the conda environment. Run this to quit the QIIME2 environment and return to the base environment.
    
        ```bash
        conda deactivate
        ```

- **MAC OS (Apple Silicon)**

    1. Update conda to the latest version.
    
        ```bash
        conda update conda
        ```
    
    2. Create a new conda environment for QIIME 2 and install it. You can choose whatever name you’d like for the environment, here we call it 'qiime2-amplicon-2025.4' to indicate what QIIME2 release is installed.
    
        ```bash
        CONDA_SUBDIR=osx-64 conda env create --name qiime2-amplicon-2026.1 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
        ```
    
    3. Activate the conda environment containing QIIME2. You will need to run this command every time you need to use QIIME2.
    
        ```bash
        conda activate qiime2-amplicon-2026.1
        conda config --env --set subdir osx-64
        ```
    
    4. Deactivate the conda environment. Run this to quit the QIIME2 environment and return to the base environment.
    
        ```bash
        conda deactivate
        ```



#### 1.1.4. Install optional softwares

We will use a couple of additional software tools to generate a report on raw sequence quality. If you wish, you can install and run them yourself; otherwise, you can rely on the report files that will be provoided.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://seqera.io/multiqc/)

---

### 1.2. Download data
You will learn to perform metabarcoding analysis using a small subset of the data collected by EUROSTAT during the Land Use Change Area Survey ([LUCAS](https://ec.europa.eu/eurostat/statistics-explained/index.php?title=LUCAS_-_Land_use_and_land_cover_survey)). Within this project soil samples were collected across Europe to characterize both soil physico-chemical properties and the biological communities present (bacteria, fungi, metazoan). Since the focus of this course are invertebrates we will use metabarcoding data produced targeting the 18S marker gene provided in this [paper](https://doi.org/10.1111/gcb.16871). Below you can see a map of all the LUCAS sampling points.

![LUCAS sampling points.](https://github.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/blob/main/img/LUCAS_sampling.jpg?raw=true)

You can download the subsample of LUCAS data and metadata that we will use in the analisis [here..missing.link..](missing.link)
> [!NOTE]
> STILL HAVE TO ADD LINK TO DOWNLOAD DATA AND METADATA

---






