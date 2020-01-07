#ifndef PRINTER_H_
#define PRINTER_H_


#include "PPMutil.h"
#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"
#include "PPMeos.h"
#include "PPMstatus.h"


#include <string>
#include <iostream>
//#include <vector>
//#include <math.h>






class Printer{
private:
    
    double rapidity_window;
    double transverse_square;
    
    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<EOS> eos;
    std::shared_ptr<Coordinates> coord;
    
    int run_num;
    
    int grid_nx, grid_ny, grid_neta;
    
    std::string profile_filename;
    
    std::string GenerateProfileFilename();
    
    void SurvayConfiguration();
    
    double hbc3 = 1.0;
    
public:
    Printer( int run_num_in,
            std::shared_ptr<EOS> eos_in,
            std::shared_ptr<Coordinates> coord_in,
            InitData &DATA_in,
            SCGrid &arena_in );
    ~Printer();// destructor
    void PrintProfile( int ppm_status );
    
    
    
};

#endif 
