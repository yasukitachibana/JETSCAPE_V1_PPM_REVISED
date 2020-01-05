#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <cmath>
#include <iostream>
#include <MakeUniqueHelper.h>

#include "JetScapeLogger.h"
#include "PPM.h"

using namespace Jetscape;

PPM::PPM() : FluidDynamics() {
    hydro_status = NOT_START;
    run_num = 0;
    SetId("PPM");
    VERBOSE(8);
}


PPM::~PPM() {
    VERBOSE(8);
}

void PPM::InitTask()
{
    VERBOSE(8);
    //kind of stupid ... do pointer GetHydroXML() via XML instance ...
}

void PPM::InitializeHydro(Parameter parameter_list) {
    
    VERBOSE(8);
    Welcome();
    
    DATA.whichEOS = 1;
    DATA.profileType = 4;
    
    DATA.nt = 40;
    DATA.nx = 193;
    DATA.ny = 193;
    DATA.neta = 95;
    DATA.T0 = 0.500;
    DATA.tau0 = 0.6;

    DATA.delta_tau = 0.3;
    DATA.delta_x = 0.3;
    DATA.delta_y = 0.3;
    DATA.delta_eta = 0.3;
    
    DATA.source = 2;
    
    DATA.surface_filename_head = "surface";
    
    DATA.profile_input_file = "/Users/yasukitachibana/Dropbox/Codes/JETSCAPE_V1_PPM_SRC_ADS/hydro_profile/Smooth_initial.dat";
    
    
    //setup EOS
    eos = nullptr;
    switch(DATA.whichEOS){
        case 1: eos = std::make_shared<BMW>(); break;
        case 0: eos = std::make_shared<MasslessIdeal>(); break;
        default: eos = std::make_shared<EOS>();
    }

    //setup Initial
    initial = std::unique_ptr<Initial> (new Initial( eos, DATA, arena ));

}

void PPM::SetInitialProfile(){
    
    if( DATA.profileType == 0 && pre_eq_ptr != nullptr ) {
        double dx = ini->GetXStep();
        double dz = ini->GetZStep();
        double z_max  = ini->GetZMax();
        int nz = ini->GetZSize();
        initial->get_preequilibrium_vectors(dx, dz, z_max, nz,
                                            pre_eq_ptr->e_,
                                            pre_eq_ptr->utau_, pre_eq_ptr->ux_,
                                            pre_eq_ptr->uy_,   pre_eq_ptr->ueta_,
                                            pre_eq_ptr->pi00_, pre_eq_ptr->pi01_, pre_eq_ptr->pi02_,
                                            pre_eq_ptr->pi03_, pre_eq_ptr->pi11_, pre_eq_ptr->pi12_,
                                            pre_eq_ptr->pi13_, pre_eq_ptr->pi22_, pre_eq_ptr->pi23_,
                                            pre_eq_ptr->pi33_, pre_eq_ptr->bulk_Pi_);
        
    }else if( DATA.profileType == 0 ){
        JSWARN<< "<-[PPM] Please Set Initial Condition ->";
        exit(-1);
    }
    
    initial->InitArena( run_num );
    
}


void PPM::EvolveHydro() {
    VERBOSE(8);
    SetInitialProfile();
    hydro_status = INITIALIZED;
    RunHydro();
    hydro_status = FINISHED;
    run_num++;
}

void PPM::RunHydro() {

    evolve = std::shared_ptr<Evolve> (new Evolve( run_num, eos, DATA, arena ));
    
    evolve->EvolveIt();
    
}

void PPM::GetHydroInfo(real t, real x, real y, real z,
                       //                           FluidCellInfo* fluid_cell_info_ptr) {
                       std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr){
    // create the unique FluidCellInfo here
    fluid_cell_info_ptr=make_unique<FluidCellInfo>();
    
    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities
    
    if (hydro_status == FINISHED)
    {
        JSWARN<< "<-[PPM] No hydro info transfer ... ->";
        exit(-1);
    }
    else
    {
        JSWARN<< "<-[PPM] Hydro not run yet ... ->";
        exit(-1);
    }
}


void PPM::Welcome(){
    
    JSINFO << "                                             ";
    JSINFO << "=============================================";
    JSINFO << "                                             ";
    JSINFO << "      -------------------------------------  ";
    JSINFO << "     / ⚠ CAUTION ⚠                        /  ";
    JSINFO << "    / This JETSCAPE is HACKED with       /   ";
    JSINFO << "   / Numerical Hydrodynamics [PPM]      /    ";
    JSINFO << "  / w/ Discritized Christoffel Symbols /     ";
    JSINFO << "  -------------------------------------      ";
    JSINFO << "        O                                    ";
    JSINFO << "           o                                 ";
    JSINFO << "             o   ▕▔▔▔▔▔▔▔▔▔▔▔╲               ";
    JSINFO << "                 ▕╮╭┻┻╮╭┻┻╮╭▕╮╲              ";
    JSINFO << "                 ▕╯┃╭╮┃┃╭╮┃╰▕╯╭▏             ";
    JSINFO << "                 ▕╭┻┻┻┛┗┻┻┛ ╰▏ ▏             ";
    JSINFO << "                 ▕╰━━━┓┈┈┈╭╮▕╭╮▏             ";
    JSINFO << "                 ▕╭╮╰┳┳┳┳╯╰╯▕╰╯▏             ";
    JSINFO << "                 ▕╰╯┈┗┛┗┛┈╭╮▕╮┈▏             ";
    JSINFO << "                                             ";
    JSINFO << "           [© Yasuki Tachibana]              ";
    JSINFO << "=============================================";
    JSINFO << "                                             ";
    
}
