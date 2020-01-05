#include "PPMprinter.h"
#include "JetScapeLogger.h"

using namespace Jetscape;
using namespace std;

Printer::Printer( int run_num_in,
                 std::shared_ptr<EOS> eos_in,
                 std::shared_ptr<Coordinates> coord_in,
                 InitData &DATA_in,
                 SCGrid &arena_in
                 ):DATA(DATA_in), arena(arena_in){
    
    JSINFO << "<-[PPM] Creating Printer ->";
    
    run_num = run_num_in;
    coord = coord_in;
    
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
    
    
}

Printer::~Printer(){
}// destructor
