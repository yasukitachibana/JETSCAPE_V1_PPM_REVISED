#ifndef LIQUEFIER_H_
#define LIQUEFIER_H_

//#include "PPMdata.h"
//#include "PPMgrid.h"
//#include "PPMutil.h"
//#include "PPMeos.h"
//#include "PPMhydroinfo.h"

#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"
#include "PPMfluidVal.h"
#include "PPMpartonCloud.h"

//#include "PPMutil.h"
//#include "PPMdata.h"
//#include "PPMgrid.h"



#include <vector>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>


class Evolve;


class Liquefier{
private:
    
    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<Coordinates> coord;
    std::shared_ptr<FluidValuables> fval;
    std::unique_ptr<PartonCloud> p_cloud;
    
    int run_num;
    
    std::string source_input_file;

    double t_damp;
    double tau_relax;
    double d_diff;
    
    double dr_sub;
    double subgrid_disc;
    
    double gamma_relax;
    double c_diff;

    double partonE, partonPx, partonPy, partonPz;
    double sourceE, sourcePx, sourcePy, sourcePz;

    double total_weight;
    
    double CausalDiffusionGreenFunctionSmooth(double delta_t, double delta_r);
    double CausalDiffusionGreenFunctionDeltaValue(double delta_t, double delta_r);
    void ConservationCheck();

    void AddSourceAtAPoint(PartonCloud::Droplet drop,
                           double weight,
                           double r, double theta, double phi);

public:
    Liquefier( int run_num_in,
              std::shared_ptr<Coordinates> coord_in,
              std::shared_ptr<FluidValuables> fval_in,
              InitData &DATA_in,
              SCGrid &arena_in);
    
    ~Liquefier();// destructor
    void AddSourceCart();
};

#endif 
