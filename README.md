# Numerical Code to Solve Relativistic Ideal Hydrodynamic Equation by Piecewise Parabolic Method (PPM) in $\tau$-$\eta_{\rm s}$ Coordinates with Discretized Christoffel Symbols 

## Preparation

To dowonload data file for Equatin of State, 
move to ```PPM``` and run ```getEOSforPPM.sh```:

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
you can find example executable file ```PPM_ADS``` for ideal hydro simulation and xml file ```source_ads_test.xml``` for its setting.
Those are generated from ```PPM_ADS.cc``` and ```source_ads_test.xml``` 
in ```example``` directory. 
You can run the code for one event by typing,

```bash
	./PPM_ADS source_ads_test.xml 1
```

Once finishing the run, you can find output files, 
```ppm_profile_run_0.txt``` for evolving hydro profile in binary and ```ppm_hypersurface_run_0.txt``` for freezeout hypersurface information in txt. 

## NOTES
Although the hydro codes are embedded in JETSCAPE, currently there is no communication between jet evolution and hydro evolution.
Also it should be noted that the codes to generate source term belongs to the hydro module ```PPM``` and are NOT considered as independent module here.

## Source codes for hydro
You can find ```PPM.cc``` and ```PPM.h``` in ```src/hydro``` and src codes associating them in ```src/hydro/PPM```. Those are the codes for the ideal hydro used in ```PPM_ADS.cc```
Blief explanations for some of them are below:

1. ```PPM/PPMprinter.cc```
The code to generate the output file storing the evolving hydro profile.

1. ```PPM/PPMfreezeout.c```
The code to generate the output file storing the freezeout surface information.

1. ```PPM/PPMpartonAdS.cc```
The code to load information of partons depositing energy and momentum into the medium fluid.

1. ```PPM/PPMpartonAdS.cc```
The code to give gaussian profile for the source generated for partons loaded in "PPM/PPMpartonAdS.cc".

1. ```PPM/PPMinitial.cc```
The code to generate initial condition of hydro profile. 

## XML file and Settings
Since the structure is stolen from JETSCAPE, it is the same to JETCAPE codes essentially. 

### XML File: ```source_ads_test.xml```

#### 1. ```<name>PPM</name>```
```<taus>0.6</taus>```: 



        
        <store_info>0</store_info>
        <s_factor>57.00</s_factor>
        <whichEOS>1</whichEOS><!-- (1:Lattice, 0:Ideal) -->
        <EOSfiles>../src/hydro/PPM/EOS</EOSfiles>
        
        <profileType>4</profileType><!-- (0:Input[tau-eta] 1:Bjorken[tau-eta], 2:DynamicalBrick, 3:3D Gaussian[t-z] 4:Transverse Input (from Dani) [tau-eta])-->
        <addCell>1</addCell><!-- (0:Use Cell Number in IS Module, 1:Change Cell Number in this PPM Module) -->
        <T0>0.5</T0><!-- (for Bjorken, Brick, and Gaussian) -->
        <nt>69</nt><!-- (Maximum Time Step) -->
        <dtau>0.3</dtau>
        <nx>193</nx>
        <ny>193</ny>
        <dx>0.3</dx>
        <neta>95</neta>
        <deta>0.3</deta>
        <source>2</source><!-- 0:No source, 1:Causal Diff ( Abail. Only for Cartesian Now), 2:Gaussian -->
        <profileInput>../hydro_profile/Smooth_initial.dat</profileInput><!-- Filename for Transverse Input (from Dani) -->
        <initProfileLong>1</initProfileLong><!-- 0: Bjorken, 1:Flat+Gauss Longitudinal Profile for Transverse Input (from Dani) -->
        
        <writeOutput>1</writeOutput>
        <profileOutput>ppm_profile</profileOutput>
        
        <freezeout>1</freezeout><!-- (0:No Freezeout, 1:Iso thermal, 2, Isochronous) -->
        <T_freezeout>0.14</T_freezeout>
        <surface_name>ppm_hypersurface</surface_name>
        
        <rapidity_window>3.5</rapidity_window>
        <transverse_square>12</transverse_square>
        
        <PPMsourceGauss>
          <name>PPMsourceGauss</name>
          <sourceInput>../source_list</sourceInput>
          <tau_thermal>0.0</tau_thermal>
          <sigma_trans>0.4</sigma_trans>
          <sigma_long>0.4</sigma_long>
        </PPMsourceGauss>
        
        <PPMLiquefier>
          <name>PPMLiquefier</name>
          <d_diff>0.8</d_diff>
          <t_damp>1.2</t_damp>
          <tau_relax>1.0</tau_relax>
          <sourceInput>/Users/yasukitachibana/Dropbox/Codes/JETSCAPE_V1_PPM_REVISED/build/source_ev1.txt</sourceInput>
        </PPMLiquefier>
        
        
        
    </PPM>







