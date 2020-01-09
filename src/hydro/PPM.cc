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
    tinyxml2::XMLElement *ppm=GetHydroXML()->FirstChildElement("PPM");
    
    
    if (ppm) {
        string s = ppm->FirstChildElement( "name" )->GetText();
        JSDEBUG << s << " to be initilizied ...";
        Welcome();
        
 
        double taus;
        ppm->FirstChildElement("taus")->QueryDoubleText(&taus);
        DATA.tau0 = taus;
        
        int info_memory;
        ppm->FirstChildElement("store_info")->QueryIntText(&info_memory);
        DATA.store_hydro_info_in_memory = info_memory;
        
        int weos;
        ppm->FirstChildElement("whichEOS")->QueryIntText(&weos);
        DATA.whichEOS = weos;
        if(DATA.whichEOS == 1){
            string eos_files;
            eos_files = ppm->FirstChildElement("EOSfiles")->GetText();
            DATA.eos_files = eos_files;
        }

        
        int profile_type;
        ppm->FirstChildElement("profileType")->QueryIntText(&profile_type);
        DATA.profileType = profile_type;
        
        if( DATA.profileType == 0 ){
            double s_factor;
            ppm->FirstChildElement("s_factor")->QueryDoubleText(&s_factor);
            DATA.sFactor = s_factor;
        }else if( DATA.profileType ==4 ){
            string input_profile;
            int init_profile_long;
            input_profile = ppm->FirstChildElement("profileInput")->GetText();
            ppm->FirstChildElement("initProfileLong")
            ->QueryIntText(&init_profile_long);
            DATA.profile_input_file = input_profile;
            DATA.init_profile_long = init_profile_long;
        }else if( DATA.profileType==1 ||
                  DATA.profileType==2 ||
                  DATA.profileType==3 ){
            double T0;
            ppm->FirstChildElement("T0")->QueryDoubleText(&T0);
            DATA.T0 = T0;
        }
        

        int add_cell;
        ppm->FirstChildElement("addCell")->QueryIntText(&add_cell);
        DATA.addCell = add_cell;
        
        int nt;
        ppm->FirstChildElement("nt")->QueryIntText(&nt);
        DATA.nt = nt;

        double dtau;
        ppm->FirstChildElement("dtau")->QueryDoubleText(&dtau);
        DATA.delta_tau = dtau;

        
        if( DATA.addCell==1 || DATA.profileType != 0 ){
            
            int nt,nx,ny,neta;
            ppm->FirstChildElement("nt")->QueryIntText(&nt);
            ppm->FirstChildElement("nx")->QueryIntText(&nx);
            ppm->FirstChildElement("ny")->QueryIntText(&ny);
            ppm->FirstChildElement("neta")->QueryIntText(&neta);
            DATA.nx = nx;
            DATA.ny = ny;
            DATA.neta = neta;
            
            double dx,deta;
            ppm->FirstChildElement("dx")->QueryDoubleText(&dx);
            ppm->FirstChildElement("deta")->QueryDoubleText(&deta);
            DATA.delta_x = dx;
            DATA.delta_y = dx;
            DATA.delta_eta = deta;
            
            DATA.x_size = DATA.delta_x*(DATA.nx - 1);
            DATA.y_size = DATA.delta_y*(DATA.ny - 1);
            DATA.eta_size = DATA.delta_eta*(DATA.neta - 1);
            
        }

        int source;
        ppm->FirstChildElement("source")->QueryIntText(&source);
        DATA.source = source;

        
        int write_output;
        ppm->FirstChildElement("writeOutput")->QueryIntText(&write_output);
        DATA.write_output = write_output;
        if(write_output == 1){
            string profile_output;
            profile_output = ppm->FirstChildElement("profileOutput")->GetText();
            DATA.profile_output = profile_output;
        }
        
        int freezeout;
        ppm->FirstChildElement("freezeout")->QueryIntText(&freezeout);
        DATA.fo_type = freezeout;
        if(freezeout != 0){
            double t_fo;
            string fo_surface;
            ppm->FirstChildElement("T_freezeout")->QueryDoubleText(&t_fo);
            fo_surface = ppm->FirstChildElement("surface_name")->GetText();
            DATA.temp_fo = t_fo;
            DATA.fo_surface = fo_surface;
            int surf_check;
            ppm->FirstChildElement("surface_check")->QueryIntText(&surf_check);
            DATA.surface_check = surf_check;
        }
        
        double t_sq, rap_wid;
        ppm->FirstChildElement("rapidity_window")->QueryDoubleText(&rap_wid);
        ppm->FirstChildElement("transverse_square")->QueryDoubleText(&t_sq);
        DATA.rapidity_window = rap_wid;
        DATA.transverse_square = t_sq;
        
    } else {
        JSWARN << " : PPM not properly initialized in XML file ...";
        exit(-1);
    }
    
    //setup EOS
    eos = nullptr;
    switch(DATA.whichEOS){
        case 1: eos = std::make_shared<BMW>(DATA); break;
        case 0: eos = std::make_shared<MasslessIdeal>(); break;
        default: eos = std::make_shared<EOS>();
    }

    //setup Initial
    initial = nullptr;
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
    JSINFO << "  / w/ Discretized Christoffel Symbols /     ";
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
