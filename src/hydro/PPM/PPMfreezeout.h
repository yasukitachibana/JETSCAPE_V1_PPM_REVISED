#ifndef FREEZEOUT_H_
#define FREEZEOUT_H_


#include "PPMutil.h"
#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"


#include <string>
#include <iostream>
//#include <vector>
//#include <math.h>






class Freezeout{
private:
    
    double temp_fo;

    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<Coordinates> coord;

    int run_num;
    
    int grid_nx, grid_ny, grid_neta;
    
    std::string surface_filename;

    std::string GenerateSurfaceFilename();


public:
    Freezeout( int run_num_in,
              std::shared_ptr<Coordinates> coord_in,
              InitData &DATA_in,
              SCGrid &arena_in );
    ~Freezeout();// destructor
    

    
    
    
    
    
//
//    int FindFreezeoutSurface(Grid ***arena);


    
};

#endif 
