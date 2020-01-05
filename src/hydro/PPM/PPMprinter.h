#ifndef PRINTER_H_
#define PRINTER_H_


#include "PPMutil.h"
#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"
#include "PPMeos.h"


#include <string>
#include <iostream>
//#include <vector>
//#include <math.h>






class Printer{
private:
    
    double temp_fo;
    
    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<EOS> eos;
    std::shared_ptr<Coordinates> coord;
    
    int run_num;
    
    int grid_nx, grid_ny, grid_neta;
    
    
public:
    Printer( int run_num_in,
            std::shared_ptr<EOS> eos_in,
            std::shared_ptr<Coordinates> coord_in,
            InitData &DATA_in,
            SCGrid &arena_in );
    ~Printer();// destructor
    
    
    
};

#endif 
