# Installation Guide

This guide will help you set up the environment and install the necessary tools for using the MonViso Reloaded package.

## Requirements

Before proceeding with the installation, ensure you have the following prerequisites:

- **Python 3.9:** The core language used for the project.
- **Cobalt:** A multiple sequence alignment tool. Download from [NCBI Cobalt](https://www.ncbi.nlm.nih.gov/tools/cobalt/cobalt.cgi?CMD=Doc).
- **HMMER:** Software for searching sequence databases for sequence homologs. Available at [HMMER](http://hmmer.org/).
- **Modeller:** Used for homology or comparative modeling of protein three-dimensional structures. Download from [Modeller](https://salilab.org/modeller/download_installation.html).

## Installation Steps

Follow these steps to install the necessary tools and the MonViso Reloaded package:

### 1. Install Dependencies

- **HMMER:** Follow the installation instructions provided on the HMMER website.
- **COBALT:** Refer to the Cobalt download page for detailed installation instructions.
- **MODELLER:** Visit the Modeller website for installation guides tailored to your operating system.

### 2. Check Python Version

Ensure you have Python 3.9 installed by running the following command in your terminal:

```bash
python --version
```

### 3. Clone the GitHub Repository
Clone the MonViso Reloaded repository to your local machine:

```bash
git clone https://github.com/alisamalb/monviso_reloaded
cd monviso_reloaded
```
### 4. Run the Setup
Install the package by executing:

```bash
python setup.py install
```

## Verifying the installation
After completing the installation steps, verify that the MonViso Reloaded package and all dependencies are correctly installed by running a quick test command or script provided by the package documentation.