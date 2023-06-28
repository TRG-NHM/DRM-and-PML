# DRM and PML

This repo is the implementation of the Domain Reduction Method (DRM) and Perfectly Matched Layer (PML) for Abaqus in Python. Although it is possible to utilize parts of the function and use it with another Finite Element Analysis (FEA) software, the development of this repo focuses on integration with Abaqus.

## Installation and Requirements

There is no need to "install" this repo. Simply clone this repo via git commands or download it from this webpage as a zip file and extract it on your machine.

However, to use this repo, you would need [Python](https://www.python.org/) 3 (>=3.10) and other required libraries. All required libraries are listed in the "requirements.txt" file, and you can easily install them with `pip`. 

    python3 -m pip install -r requirements.txt

Optionally, you would need Open MPI if you want to take advantage of parallel computing for specific functionalities. Note that you may also need to install GNU Scientific Library (GSL) to make Open MPI work properly. For UNIX-like systems (e.g., macOS), the easiest way to install Open MPI and GSL is by using [Homebrew](https://brew.sh/). Assuming Homebrew is installed on your machine, all you need to do is typing the following commands in the terminal one by one.

    brew install open-mpi
    brew install gsl

For Windows systems, it might be worth taking a look at [Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi).

## Basic Usage

The primary file you need to modify is the `__main__.py` under the `drm` folder. There are two examples (written as two functions) listed in this file for your reference. After necessary modifications that fit your need, you can run it by the following command (assuming you have already changed your current working directory to the root folder of this repo):

    python3 drm -s stepNum

`stepNum` should be changed to the step number (1, 2, 3, or 4) you want to run. The explanation for each step is included in the comment section in each function.

## Acknowledgement

The DRM and PML implementations are initially developed by Wenyang Zhang. His original implementation in MATLAB can be found [here](https://github.com/etacir/TRG-Regional-Seismic-Workflow).