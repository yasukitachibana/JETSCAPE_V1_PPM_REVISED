#ifndef SOURCEGAUSS_H_
#define SOURCEGAUSS_H_


#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMcoord.h"
#include "PPMutil.h"
#include "PPMfluidVal.h"
#include "PPMpartonAdS.h"


#include <vector>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>


class SourceGauss{
private:

    InitData &DATA;
    SCGrid &arena;
    
    std::shared_ptr<Coordinates> coord;
    std::shared_ptr<FluidValuables> fval;
    std::unique_ptr<PartonAdS> p_ads;
    
    int run_num;
    
    std::string source_input_file;
    double tau_th;
    double sigma_l;
    double sigma_t;

    double partonE, partonPx, partonPy, partonPz;
    double sourceE, sourcePx, sourcePy, sourcePz;

    double correction;

    double sigma_cut;
    double dr_sub;
    double deta_sub;
    double subgrid_disc;

    double r_max;
    int n_r;
    double dr;

    double eta_max;
    int n_eta;
    double deta;

    double Correction();
    
    void AddSourceAtAPoint(PartonAdS::Droplet drop,
                           double weight,
                           double r, double phi, double eta);

public:
    SourceGauss( int run_num_in,
    std::shared_ptr<Coordinates> coord_in,
    std::shared_ptr<FluidValuables> fval_in,
    InitData &DATA_in,
    SCGrid &arena_in );//, InitData &DATA_in, SCGrid &arena_in );
    ~SourceGauss();// destructor
    void AddSourceTauEta();
    
};

#endif 
