#ifndef INITIAL_H
#define INITIAL_H

#include <vector>
#include <cmath>
#include <float.h>
#include <array>

#include "PPMdata.h"
#include "PPMeos.h"
#include "PPMgrid.h"
#include "PPMutil.h"

class Initial
{
private:
    


    InitData &DATA;
    SCGrid &arena;
    std::shared_ptr<EOS> eos;
    
    void GenerateArena();
    
    void AddEdges();
    void SetBoundary(int extra_nx, int extra_ny, int extra_neta );
    void CopyThermo( std::array<int, 3> i_copy, std::array<int, 3> i_og );
    double GetX(int ix);
    double GetY(int iy);
    double GetEta(int ieta);

    int extra_nx;
    int extra_ny;
    int extra_neta;
    void Reshape();
    void TrentoInit();
    void InitProfileFromPreeq();
    void clean_up_arrays();
    
    void InitProfileBjorken();
    void InitProfile3DGaussian();
    
    static const int n_x_input = 201; // Dani's 2D profile table
    const double dx_input = 0.1; // Dani's 2D profile table
    const double size_input = 10.0; // Dani's 2D profile table
    void InitProfile2DFromFile( int run_num );
    void LoadInitProfile2D
    (std::array<std::array<double, n_x_input>, n_x_input> &e_trans );
    double ProfileEtaFlatGauss(double eta,
                               double eta_flat, double sigma_eta,
                               double e_cm );
    int GetIx2DInput(double x);
    int GetIy2DInput(double y);
    
    double initial_dx, initial_dz, initial_z_max, initial_nz;

    std::vector<double> initial_energy_density;
    std::vector<double> initial_u_tau;
    std::vector<double> initial_u_x;
    std::vector<double> initial_u_y;
    std::vector<double> initial_u_eta;
    std::vector<double> initial_pi_00;
    std::vector<double> initial_pi_01;
    std::vector<double> initial_pi_02;
    std::vector<double> initial_pi_03;
    std::vector<double> initial_pi_11;
    std::vector<double> initial_pi_12;
    std::vector<double> initial_pi_13;
    std::vector<double> initial_pi_22;
    std::vector<double> initial_pi_23;
    std::vector<double> initial_pi_33;
    std::vector<double> initial_bulk_pi;
    
    
public:
    Initial(std::shared_ptr<EOS> eos_in, InitData &DATA_in, SCGrid &arena_in );
    ~Initial();  //destructor

    void InitArena( int run_num );
    void get_preequilibrium_vectors( const double dx, const double dz,
                                    const double z_max, const int nz,
                                    std::vector<double> e_in,
                                    std::vector<double> u_tau_in,
                                    std::vector<double> u_x_in,
                                    std::vector<double> u_y_in,
                                    std::vector<double> u_eta_in,
                                    std::vector<double> pi_00_in,
                                    std::vector<double> pi_01_in,
                                    std::vector<double> pi_02_in,
                                    std::vector<double> pi_03_in,
                                    std::vector<double> pi_11_in,
                                    std::vector<double> pi_12_in,
                                    std::vector<double> pi_13_in,
                                    std::vector<double> pi_22_in,
                                    std::vector<double> pi_23_in,
                                    std::vector<double> pi_33_in,
                                    std::vector<double> Bulk_pi_in);

};

#endif
