#include "PPMfreezeout.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace std;

Freezeout::Freezeout( int run_num_in,
                     std::shared_ptr<Coordinates> coord_in,
                     InitData &DATA_in,
                     SCGrid &arena_in
                     ):DATA(DATA_in), arena(arena_in){
    
    JSINFO << "<-[PPM] Creating Freezeout ->";
    
    run_num = run_num_in;
    coord = coord_in;
    temp_fo = DATA.temp_fo;
    
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
    
    if(DATA.fo_type == 0){
        JSINFO << "<-[PPM] Freezeout Temperature: " << temp_fo <<" GeV ->";
    }else{
        JSINFO << "<-[PPM] Isochronous Freezeout ->";
    }
    
    surface_filename = GenerateSurfaceFilename();
    
    JSINFO
    << "<-[PPM] Surface Info File Name: '"
    << surface_filename
    << "' ->";
    
}

Freezeout::~Freezeout(){
}// destructor

std::string Freezeout::GenerateSurfaceFilename(){
    return DATA.surface_filename_head + "_run_" + to_string(run_num) + ".txt";
}





//
//
//int Freezeout::FindFreezeoutSurface(Grid ***arena){
//    JSINFO << "Finding Hyper Surface --Not Implemented Yet--" ;
//    double temp_max = 0.0;
////    return 1;
//    return 0;
//}
