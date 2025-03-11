# A multi-dimensional view of a unified model for TDEs

This repository contains the parameter files and model grids to re-create
the results from the paper. It also includes the synthetic spectra and tables
of the physical parameters for all of the models shown.

- The model grids and related Python scripts are stored in the `make-grids`
  directory
- The data used in the paper are kept in the `models_1d` and `models_2d`
  directories. Note that some models include a `data` directory, which is
  atomic data not part of the default Python atomic data.
- To create the figures, use the scripts provided in `make-figures`. Make sure
  you do a recursive clone to clone one of the dependencies (which is a
  submodule in this repository).

All models were run using Python V88, commit [#4277688](https://github.com/agnwinds/python/tree/4277688729ccadc43eaf49bc7f05142a6b3c355c).

Plots were created using [PySi](https://github.com/saultyevil/pysi).

## Contact

If anything is unclear or broken, please contact me. My contact details can be
found on my GitHub profile page.

