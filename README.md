# croco-eddy
Notebook to process eddy simulations performed with CROCO

This repository lists all steps to prepare, run and analyze idealized eddy simulations with CROCO, as computed in
Barboni et al (2024, submitted to JAMES) : *How atmospheric forcing frequency, horizontal and vertical grid resolutions impact mesoscale eddy evolution in a numerical model*
http://dx.doi.org/10.13140/RG.2.2.30074.47047

## Preprocessing with croco_tools

In preprocess folder, run the matlab script `make_vortex.m` to create an initialized mesoscale vortex , the corresponding batch to do this on Datarmor is available.
Eddy and background parameters to adjust your eddy should be changed here. Background is set as in Eq.5 in Barboni et al (2024), and eddy equation as in Eq. 6.

## Run the CROCO model

In `croco_run/` floder, set up your configuration in `cppdefs.h` and grid/parallelization parameters in `param.h`.
The modified test case used in Barboni et al (2024) is under `VORTEX_LMD` key.

### Simulations with bulk fluxes air-sea interactions (COARE parametrization)

Keys `BULK_FLUX`, `BULK_LMD` and `BULK_SM_UPDATE` should be DEFINED , keys `ANA_STFLUX`,`ANA_SMFLUX`,`IN_REAL_WIND`,`IN_REAL_STFLUX` and `IN_REALSRFLUX` should NOT be defined

Current feedback option is available with key `CFB_STRESS`. This configuration was not tested in Barboni et al (2024), but comparison is shown in Fig.5.4 and 5.5 in Barboni, A. (2023), [PhD manuscript] Surface and subsurface evolution of mesoscale eddies under atmospheric forcing: case study in the
Mediterranean Sea
Thermal current feedback (ocean current retroaction on the thermal flux) is available (but not tested) with key `CFB_STRESS`

### Simulations with direct air-sea flux forcing (no COARE parametrization)

Keys  `ANA_STFLUX`,`ANA_SMFLUX`,`IN_REAL_WIND`,`IN_REAL_STFLUX` and `IN_REALSRFLUX` should be DEFINED , keys  `BULK_FLUX`, `BULK_LMD` and `BULK_SM_UPDATE` should NOT be defined
Current feedback option is available (but not tested)  with key `IN_WIND_CFB`

### Compilation and run
Set up in `jobcomp` the `SOURCE` variable to point to your CROCO code repository
Compile (`batch_compil` to do this on Datarmor)
A `Compile` folder should appear with raw Fortran routines (.F) and precompiled ones (_.f), check that your changes and key choices were taken into account in `_.f` files.

Set up experiment parameters (timestep, recorded and averaged variables, output path, etc.) in `croco.in.QWA`. Best practice is `expname_his.nc` for history variables and `expname_avg.nc` for averaged ones.
Execute the code with `batch_exe_croco`

## Eddy tracking postprocessing

Output files are usually very heavy. First slice the model output to retain only the *averaged* sea surface height :
`ncks -v zeta /output-path/expname_avg.nc /some-work-space-path/expname_ssh.nc`

Compute geostrophic speed using `make_geo_file_zeta.m` (you can launch it with `batch_geos_matlab`). zonal and meridional speed are computed in `expname_ssu.nc`and `expname_ssv.nc`

Then performed AMEDA eddy tracking in folder `ameda_postproc/`. User parameters are set in `keys_sources_WindEddy_.m`. AMEDA code is as in https://github.com/briaclevu/AMEDA
Some modifications (mostly grid normalization) can be tracked typing `grep -r WindEddy *`
- source should be `WindEddy`, config is `''`
- `runname` is experiment folder name, `expname` is output file name without `_ssh.nc`
- change nc_u, nc_v, nc_ssh file names accordingly, same for variable names
- WARNING Coriolis parameter `f_0` and time step `dps` should be manually set. `dps` is the time between 2 averaging steps
- `deg` factor can be increased to reduce computation time (typically `deg=2` for 1km grid size)
- `borders` set the border replication to mimic periodic borders (typically `borders=40` for `deg=2` and 1km grid size)

Launch AMEDA running `MAIN_AMEDA_WindEddy.m` (no modification needed there).

## Experiment analysis

First use the notebook `AMEDA-Vort-SST-MLD.ipynb` to compute 2 images
Then use `AMEDA-N2-MLD-vert.ipynb` to extract relevant analysis and store it in a light .npy file per experiment
Use `Eddy_timeseries_deltaT-Q-Mix.ipynb` to gather all .npy file and produce the comparison figures (Fig. 6 and 9 in Barboni et al, 2024)
Near-inertail waves spectral investigation is done in `AMEDA_Buoyancy_flux_NIIW.ipynb` (cf Fig.10 in Barboni et al, 2024)







