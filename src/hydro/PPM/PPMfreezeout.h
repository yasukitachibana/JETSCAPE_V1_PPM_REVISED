#ifndef FREEZEOUT_H_
#define FREEZEOUT_H_


#include "PPMutil.h"
#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"
#include "PPMstatus.h"
#include "PPMfluidVal.h"
#include "PPMeos.h"
#include "PPMsurfaceCheck.h"

#include <string>
#include <iostream>
#include <cfloat>
//#include <vector>
//#include <math.h>






class Freezeout{
private:
    
    double temp_fo;
    double rapidity_window;
    
    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<Coordinates> coord;
    std::shared_ptr<FluidValuables> fval;
    std::shared_ptr<EOS> eos;
    
    int run_num;
    
    int grid_nx, grid_ny, grid_neta;
    
    std::string surface_filename;
    
    std::string GenerateSurfaceFilename();
    
    std::ofstream ofs_freezeout;
    void OpenSurfaceFile();
    
    void IsochronousFreezeout( int threshold );
    int FullFreezeout();
    void BulkFreezeout(const std::array<int, 3> &i_cell,
                       const std::array<double, 3> &x_cell,
                       const std::array<double, 4> &dsigma);
    void SurfaceFreezeout(const std::array<int, 3> &i_cell,
                          const std::array<double, 3> &x_cell,
                          const std::array<double, 4> &dsigma );
    std::array<double, 4> GetDsigma( int time_shift );
    
    std::array<std::array<int, 3>, 3> point_next{{{1,0,0},
                                                  {0,1,0},
                                                  {0,0,1}}};
    
    double GetE( double U, double u0, double ux );
    
public:
    Freezeout( int run_num_in,
              std::shared_ptr<Coordinates> coord_in,
              std::shared_ptr<FluidValuables> fval_in,
              std::shared_ptr<EOS> eos_in,
              InitData &DATA_in,
              SCGrid &arena_in );
    ~Freezeout();// destructor
    
    int FindFreezeoutSurface(int ppm_status);
    
    
    
};

#endif 
