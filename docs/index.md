# Welcome to MoNvIso Documentation

MoNvIso is a comprehensive software tool designed for the analysis and modeling of protein isoforms. It automates the process of identifying canonical and additional isoforms, assessing their modeling propensity, mapping mutations accurately, and building structural models of proteins. By leveraging data from the Uniprot database, MoNvIso facilitates a deeper understanding of protein function and variation.

## Overview

The MoNvIso workflow is streamlined into three primary steps, making it both powerful and user-friendly for researchers and bioinformaticians alike. This workflow ensures a meticulous examination of gene names, isoforms, mutations, and the subsequent construction of structural models when experimental structures are not available.

### Step 1: Gene and Isoform Verification

MoNvIso begins by verifying gene names provided in the user's input file, identifying both canonical and additional isoforms directly from the Uniprot database. The input file also includes a list of mutations of interest, setting the stage for detailed analysis in subsequent steps.

### Step 2: Modelling Propensity and Mutation Mapping

This step focuses on assessing the suitability of available isoforms for modeling and how mutations map onto these isoforms. MoNvIso supervises the availability of templates and ensures mutations are associated with the appropriate isoforms. It identifies and highlights any discrepancies where mutations cannot be mapped, ensuring accuracy and reliability in the selection process.

The selection procedure, known as the Selection Function, integrates two critical components:

 - Structural Function: Assesses the modelling propensity.
 - Mutation Function: Evaluates the mapping of mutations on the isoforms.

Both functions are weighted (default: w1 = w2 = 10), which can be adjusted by the user to tailor the analysis.

### Step 3: Structural Model Building

In the final step, MoNvIso builds structural models of the identified proteins, including both wild-type (WT) forms and their variants, selected based on the analysis conducted in Step 2. This step is only initiated if experimental structures are not already available, ensuring that only necessary modeling work is performed.
Getting Started

Dive into the world of protein isoform analysis with MoNvIso by exploring the following sections of our documentation:

- [Installation Guide](installation.md): Learn how to install MoNvIso on your system.
- [Quick Start Tutorial](tutorial.md): Get up and running with your first analysis.
- [Understanding the Workflow:](workflow.md) A deeper dive into each step of the MoNvIso workflow.
- [FAQs](FAQs.md): Answers to common questions about MoNvIso.

MoNvIso is here to enhance your research with advanced protein isoform analysis and modeling. Whether you're investigating mutations, studying protein functions, or exploring structural variations, MoNvIso provides the tools you need to achieve your objectives with precision and ease.
