# Example script for PyDePSI

This folder contains scripts to run PyDePSI on SLURM infrastructures. One needs to modify `exec_python.slurm` to specify 1) a conda environment with PyDePSI and 2) the desired Python script. The command can be used to execute the specified script:

```bash
sbatch exec_python.slurm
```

The submitted SLURM job will execute the Python script. A Dask cluster will be created, with the submitted job as the scheduler. More SLURM jobs will be submitted from the first job to start the Dask workers. The Dask cluster will be used to parallelize the computation of the Python script.

## Prerequisites

It is assumed that the user has a working conda environment with PyDePSI installed. This conda environment is required in the `exec_python.slurm` script. 

Some HPC system (for example SURF Spider system, see [this documentation](https://doc.spider.surfsara.nl/en/latest/Pages/software_on_spider.html)) may require the user to execute the conda environment inside a container. One can for example use the [LUMI container wrapper](https://docs.lumi-supercomputer.eu/software/installing/container-wrapper/) as a solution.

## PyDePSI with RS-DAT

The [RS-DAT platform](https://rs-dat.github.io/RS-DAT/) developed by Netherlands eScience Center can be used to start up a Jupyter environment with Dask cluster. With RS-DAT, one can also start a Dask cluster schedular on a HPC system, then pass the schedular address to the Python script. In this way, the same Dask cluster can be re-used for multiple Python script executions. This can be useful for developing and testing PyDePSI scripts.

Please refer to the [RS-DAT JupyterDaskOnSLURM user guide](https://github.com/RS-DAT/JupyterDaskOnSLURM/blob/main/user-guide.md) for this option.