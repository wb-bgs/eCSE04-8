help([[
WMAM 3.7
========
The World Magnetic Anomaly Modelling (WMAM) application is an MPI code for mapping the natural
magnetisation of the Earth's surface rocks. This crustal field map takes the form of a spherical
harmonic model.

The name of the WMAM executable file is "mod_wmam_020"; the location of which is placed on the path.

The executable file takes four arguments, the spherical harmonic degree, the resolution degree,
the inversion scheme id (Polak-Ribiere=1, Conjugate-Gradient=2) and the damping factor (usually 5.0).

Typical values for the spherical harmonic degree and resolution degree are (200,1.0), (300,0.5),
(720,0.25), (1440,0.1) and (2000,0.05).

The input data files must reside in a folder called "Data" within the Slurm submission directory.
Similarly, output data files are written to folder called "Results".

A sample submission script can be found at "/work/y07/shared/apps/core/wmam/3.7/share/submit.ll".


Three input files are required by WMAM.

1. coef_1990_15.dat.bin
The model coefficients consistent with the 1990 terrestrial core field.

2. wdmam_geocentric.dat.bin
The magnetic anomaly data.

3. model.in.bin
The starting coefficients for the output field mode.


As well as the Slurm output/error files, three files are produced by WMAM.

1. model_No_P.out
The final field model coefficients.

2. fit_No_P.out.bin
Diagnostic info relating to the geographic points covered by the anomaly data.

3. fit_damp.out.bin
Diagnostic info relating to the nodal points of the spherical harmonic model.


For Further info about the WMAM app, see "https://github.com/wb-bgs/eCSE04-8".
Please email "wb@bgs.ac.uk" for access to this repository.


   - Installed by: M. Bareford, EPCC"
   - Date: May 2023\n"

]])

load("cpe/21.09")
prepend_path("LD_LIBRARY_PATH", os.getenv("CRAY_LD_LIBRARY_PATH"))

pushenv("SLURM_CPU_FREQ_REQ", "2250000")

local pkgName = myModuleName()
local pkgNameVer = myModuleFullName()
local pkgNameBase = pathJoin("/work/y07/shared/apps/core", pkgName)
local pkgVersionBase = pathJoin("/work/y07/shared/apps/core", pkgNameVer)

prepend_path("PATH", pathJoin(pkgVersionBase, "bin"))

family("wmam")
