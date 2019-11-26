# Efficient classical simulation of random shallow 2D quantum circuits: code

This repository contains code to reproduce numerical experiments.

---

### Installation

* **Requirements**:
 macOS, [`itensor version 2`](https://itensor.org/), [`numpy`](https://www.numpy.org), [`matplotlib`](https://matplotlib.org/), [`jupyter`](https://jupyter.org/)

* **Download** the repository:
  ```
  mkdir random-shallow-2d
  git clone https://github.com/random-shallow-2d/random-shallow-2d.git random-shallow-2d
  cd random-shallow-2d
  ```
* **Compile**:
  If your installation of iTensor is not in your home directory, you will need to edit the Makefile. In particular, this line:
  ```
  LIBRARY_DIR=$(HOME)/itensor
  ```
  must be changed to the correct location of iTensor. Now you are ready to compile. Simply run:
  ```
  make
  ```
  which should create two executables: `cluster` and `brickwork`.
  
---

### Usage

The simulations are most easily run via the interface of the Jupyter notebook `interface.ipynb`:
  ```
  jupyter notebook interface.ipynb
  ```
From the notebook, the user can run simulations, analyze data, and generate plots. Examples and further usage instructions are given in the notebook.
