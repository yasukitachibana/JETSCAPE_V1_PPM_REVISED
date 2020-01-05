#ifndef PPM_H
#define PPM_H

#include "FluidDynamics.h"
#include "PPM/PPMdata.h"
#include "PPM/PPMeos.h"
#include "PPM/PPMgrid.h"
#include "PPM/PPMinitial.h"
#include "PPM/PPMevolve.h"

//#include "PPM/PPMhydroinfo.h"

using namespace Jetscape;

class PPM: public FluidDynamics {

private:

    int run_num;
    
    InitData DATA;
    std::shared_ptr<EOS> eos;
    std::unique_ptr<Initial> initial;
    std::shared_ptr<Evolve> evolve;
    
    SCGrid arena;
    
    void SetInitialProfile();
    void RunHydro();
    void Welcome();

  
public:
    
    PPM();
    ~PPM();

    void InitializeHydro(Parameter parameter_list);
    void EvolveHydro();
    void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
		      std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);
    void InitTask();
    //virtual void Exec();

     
  };

#endif  // PPM_H
