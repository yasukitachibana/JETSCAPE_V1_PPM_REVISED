# Numerical Code to Solve Relativistic Hydrodynamic Equation with Piecewise Parabolic Method (PPM)

## Preparation

To dowonload data file for Equatin of State, 
move to "PPM" and run "getEOSforPPM.sh":

```bash
	cd YOUR_CODE_LOCATION/src/hydro/PPM
	./get_eos_for_ppm.sh
	cd YOUR_CODE_LOCATION
```

## Compile and Run

With cmake:

```bash
	mkdir build
    cd build
    cmake ..
    make -j
```

If you do not have any problem in compilation,
you can find example executable file "PPM\_ADS" for ideal hydro simulation and xml file "source\_ads\_test.xml" for its setting.
Those are generated from "PPM\_ADS.cc" and "source\_ads\_test.xml" 
in "example" directory. 
You can run the code for one event by typing,

```bash
	./PPM_ADS source_ads_test.xml 1
```

Once finishing the run, you can find output files, 
"ppm\_profile\_run\_0.txt" for evolving hydro profile in binary and "ppm\_hypersurface\_run\_0.txt" for freezeout hypersurface information in txt. 

## NOTES
Although the hydro codes are embedded in JETSCAPE, currently there is no communication between jet evolution and hydro evolution.
Also it should be noted that the codes to generate source term belongs to the hydro module "PPM" and are NOT considered as independent module here.

## Source codes for hydro
You can find "PPM.cc" and "PPM.h" in "src/hydro" and src codes associating them in "src/hydro/PPM". Those are the codes for the ideal hydro used in "PPM\_ADS.cc"
Blief explanations for some of them are below:

1. "PPM/PPMprinter.cc"
The code to generate the output file storing the evolving hydro profile.

1. "PPM/PPMfreezeout.cc"
The code to generate the output file storing the freezeout surface information.

1. "PPM/PPMpartonAdS.cc"
The code to load information of partons depositing energy and momentum into the medium fluid.

1. "PPM/PPMpartonAdS.cc"
The code to give gaussian profile for the source generated for partons loaded in "PPM/PPMpartonAdS.cc".

1. "PPM/PPMinitial.cc"
The code to generate initial condition of hydro profile. 

## XML file and Settings
Since the structure is stolen from JETSCAPE, it is the same to JETCAPE codes essentially. 

1. PPM







