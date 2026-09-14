# CNV-analysis-Bertero-lab
CNV analisis for cell kariotyping

# CNV Analysis with GenomeStudio and RStudio Server

This repository provides a Docker-based workflow for CNV analysis using:

- **GenomeStudio** for SNP-array processing
- **RStudio Server** for downstream CNV analysis and visualization

---

## 1. Prepare the experiment folder

Create a folder for your experiment.

This folder will be referred to as:

```bash
CNV_ROOT
CNV_ROOT=/path/to/your/experiment_folder
```
Note: CNV_ROOT must point to your experiment folder.

---

## 2. Add the GenomeStudio startup script

Place the GenomeStudio startup script inside the experiment folder:
- **start-genomestudio.sh for macOS / Linux
- **start-genomestudio.bat for Windows
This script is used to start the GenomeStudio Docker container.

---

## 3. Create the input folder structure
Inside the experiment folder, create:

scratch/
└── input_files/
    └── genomestudio/

The scratch folder will be mounted inside the Docker containers.
Place the files required by GenomeStudio inside:

scratch/input_files/genomestudio/
