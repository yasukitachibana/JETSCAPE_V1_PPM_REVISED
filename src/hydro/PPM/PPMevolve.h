#ifndef EVOLVE_H_
#define EVOLVE_H_

#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMutil.h"
#include "PPMeos.h"
#include "PPMcoord.h"
#include "PPMfluidVal.h"
#include "PPMliquefier.h"
#include "PPMsourceGauss.h"
#include "PPMfreezeout.h"
#include "PPMprinter.h"


#include <memory>


class Evolve{
private:
    
    int run_num;
    
    int grid_nx, grid_ny, grid_neta;
    
    std::shared_ptr<EOS> eos;
    InitData &DATA;
    SCGrid &arena;
    
    int ppm_status;

    
    std::shared_ptr<Coordinates> coord;
    std::unique_ptr<Liquefier> liquefier;
    std::unique_ptr<SourceGauss> source_gauss;
    std::unique_ptr<Freezeout> freezeout;
    std::shared_ptr<FluidValuables> fval;
    std::unique_ptr<Printer> printer;
    
    
    
    void InitialSetting();
    
    double (Evolve::*F_0)(double U, double vx, double p);
    double (Evolve::*F_i)(double U, double vx, double p);
    void (Evolve::*Flux[3])(std::array<double, 5> &flux,
                            const std::array<int, 3> &i_evo_m_1,
                            const std::array<int, 3> &i_evo,
                            double dtdx);

    double F_x(double U, double vx, double p);
    double F_tau(double U, double vx, double p);
    double F_xCart(double U, double vx, double p);
    double F_tCart(double U, double vx, double p);
    double Fcurrent(double U, double vx);
    
    void FluxTransEta(std::array<double, 5> &flux,
                      const std::array<int, 3> &i_evo_m_1,
                      const std::array<int, 3> &i_evo,
                      double dtdx);
    
    void FluxTrans(std::array<double, 5> &flux,
                   const std::array<int, 3> &i_evo_m_1,
                   const std::array<int, 3> &i_evo,
                   double dtdx);
    
    void StepDtau( int it, std::array<int, 3> &evoDir );
    void MiniStep( double strang, int direction );
    
    void SetFluidEvoDir(int direction,
                        std::array<int, 3> &s, std::array<int, 3> &ncell);
    
    void SetFluidEvoIndexSet( std::array<std::array<int, 3>, 5> &i_evo,
                             int i, int j, int k,
                             const std::array<int, 3> &s);
    
    void SetFluidEvoIndex( std::array<int, 3> &i_evo,
                          int i, int j, int k,
                          const std::array<int, 3> &s);
    
    void SetInitUmem( std::array<std::array<double, 5>, 5> &U_mem,
                     const std::array<std::array<int, 3>, 5> &i_evo );
    
    void GetU_R( std::array<double, 5> &U_R,
                const std::array<double, 5> &U,
                const std::array<double, 5> &Up1,
                const std::array<double, 5> &Um1,
                const std::array<double, 5> &Up2 );
    double Dmu(double Up1, double U, double Um1);
    
    void SetMono(std::array<double, 5> &UR,
                 std::array<double, 5> &UL,
                 std::array<double, 5> &UU );
    void SetMonoComponent(double& UR, double& UL, double& UU);
    
    void GetBorderValuables(std::array<double, 5> &U, int direction,
                            double &vx, double &c, double &p);
    
    void SetUmem(const std::array<int, 3> &i_evo,
                 std::array<std::array<double, 5>, 5> &U_mem,
                 std::array<std::array<double, 5>, 2> &U_L,
                 std::array<std::array<double, 5>, 2> &U_R,
                 std::array<double, 2> &v_L,
                 std::array<double, 2> &v_R,
                 std::array<double, 2> &c_L,
                 std::array<double, 2> &c_R );
    
    double VelocitySum( double v1, double v2 );
    
    void GetUbarR( std::array<double, 5> &UbarR,
                  const std::array<double, 5> &U,
                  const std::array<double, 5> &UL,
                  const std::array<double, 5> &UR,
                  double v_r, double c_r,
                  double v_l_p, double c_l_p,
                  double dtdx);
    
    double GetUbarRComponent(double U, double UL, double UR, double  lambda);
    
    void GetUbarL( std::array<double, 5> &UbarL,
                  const std::array<double, 5> &U,
                  const std::array<double, 5> &UL,
                  const std::array<double, 5> &UR,
                  double v_l, double c_l,
                  double v_r_m, double c_r_m,
                  double dtdx);
    
    double GetUbarLComponent(double U, double UL, double UR, double lambda);
    
    void GetFluxHlle(std::array<double, 5> &flux,
                     std::array<double, 5> &U_minus,
                     std::array<double, 5> &U_plus,
                     int direction);
    
    double GetBminus(double v_h, double c_h, double v_minus, double c_minus);
    double GetBplus(double v_h, double c_h, double v_plus, double c_plus);
    
    void GetFh(std::array<double, 5> &flux,
               const std::array<double, 5> &U_plus,
               const std::array<double, 5> &U_minus,
               double b_plus, double b_minus,
               double v_plus, double v_minus,
               double p_plus, double p_minus,
               int direction);
    
    double GetFhXComponent(double U_plus, double U_minus,
                           double b_plus, double b_minus,
                           double v_plus, double v_minus,
                           double p_plus, double p_minus);
    double GetFh0Component(double U_plus, double U_minus,
                           double b_plus, double b_minus,
                           double v_plus, double v_minus,
                           double p_plus, double p_minus);
    double GetFhCurrentComponent(double U_plus, double U_minus,
                                 double b_plus, double b_minus,
                                 double v_plus, double v_minus);

    void DirExchange( int it, std::array<int, 3> &evoDir );    
    
public:
    Evolve( int run_num, std::shared_ptr<EOS> eos_in, InitData &DATA_in, SCGrid &arena_in );
    ~Evolve();// destructor
    void EvolveIt();
    
};

#endif 
