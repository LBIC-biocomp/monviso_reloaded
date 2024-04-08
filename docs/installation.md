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

### 3. Install via pip

```bash
pip install monviso
```

## Verifying the installation
After completing the installation steps, verify that the MonViso Reloaded package and all dependencies are correctly installed by running a quick test command or script provided by the package documentation.

### Alternative - Clone the GitHub Repository
Clone the MonViso Reloaded repository to your local machine:

```bash
git clone https://github.com/alisamalb/monviso_reloaded
cd monviso_reloaded
```

Install the package by executing:

```bash
python setup.py install
```

### Alternative - Run the Docker container

1. **Download the Containerfile:** First, you need to download the Containerfile. Save the file from this URL: [Containerfile](https://raw.githubusercontent.com/alisamalb/monviso_reloaded/main/Containerfile).

2. **Insert Your Modeller License:** Open the downloaded Containerfile in a text editor and replace `XXXX` with your Modeller license key. This step is crucial for ensuring that Modeller functions correctly within the container.

3. **Build the Docker Image:** With the license key in place, you can now build your Docker image. Open a terminal and run the following command in the directory where the Containerfile is located:   
```bash
docker build -t monviso -f Containerfile .
```
This command builds a Docker image named monviso using the instructions in your Containerfile.

4. **Run the Container with a Mounted Working Directory:** If you wish to use your own working directory (e.g., /Work) inside the container, you can do so by mounting it when you start the container. Run the following command to start the container and access your working directory within it:
```bash
docker run -it -v /Work:/Monviso/Work monviso
```
Replace /Work with the path to your actual working directory on the host machine. This command mounts your host directory /Work to /Monviso/Work inside the container, allowing you to work with your files directly from within the container environment.

By following these steps, you'll have a fully functional Monviso environment set up and ready to use in Docker, complete with your required Modeller license and access to your working directory.

