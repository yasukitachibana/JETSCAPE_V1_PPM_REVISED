#ifndef SURFACECHECK_H_
#define SURFACECHECK_H_

#include "PPMutil.h"

#include <string>
#include <iostream>

class SurfaceCheck{
    
private:
    std::string filename;
    
    double total_energy;
    double total_momx;
    double total_momy;
    double total_momz;
    
    void EndMessage();
    
public:
    SurfaceCheck( std::string filename_in );
    ~SurfaceCheck();// destructor
    
    void DoCheck();
    
    
    
};

#endif 
