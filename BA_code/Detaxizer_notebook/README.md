# Detaxizer Notebook

The **Detaxizer Notebook** that I have created in the space of the Maier Lab is an instruction and execution guide to understand how the `nf-core/detaxizer` pipeline functions, [click here](https://nf-co.re/detaxizer/1.0.0/) for more information. `FASTQ` files containing metagenomic data are taken to inspect specific taxa and also filter them out if desired.

To find out about other nf-core pipelines [click here](https://nf-co.re/).

The notebook helps any individual who might benefit from step-by-step explanations to process their data through this pipeline.

## Steps to use a Jupyter Notebook

To use and run the notebook flawlessly, it is necessary to have a Jupyter environment, as well as Conda to manage dependencies.

### Jupyter Notebook

The execution of a notebook can be mainly done in two ways:

* Within VSCode: The program has an integrated interface for Jupyter Notebooks. The `Jupyter` and `Remote-SSH` extensions must be installed if you work remotely.
* Standalone Jupyter: Jupyter Notebooks can also be installed on your systems if you run locally.

### Conda

To install packages and isolate project environments, you will need Conda. For installment instructions [click here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). If you are part of the M3 group, you can already create Conda environments.

From the provided YAML file, the VSCode conda environment can be created with this command:

```
bash conda env create -f Environments/VScode.yaml
```

To activate the environment:

```
conda activate VScode
```

Other YAML files for different environments are in the `./Environments` folder. The procedure is the same.

## Output of Pipeline

The pipeline generates a `MultiQC report` to demonstrate quality control, as well as filtered files that can later be used in the **Taxprofiler Notebook**.