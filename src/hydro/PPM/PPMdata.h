#ifndef SRC_DATA_H
#define SRC_DATA_H

#include <string>

typedef struct init_data {

    int addCell;
    
    double sFactor;

    int nx;
    int ny;
    int neta;
    int nt;

    double x_size;       /*!< in fermi -x_size/2 < x < x_size/2 */
    double y_size;       /*!< in fermi, ditto */
    double eta_size;     /*!< ditto */
    double tau_size;     /*!< tau_0 < tau < tau0+tau_size  */
    
    double delta_x;
    double delta_y;
    double delta_eta;
    double delta_tau;
    
    int store_hydro_info_in_memory;

    double tau0;
    int whichEOS; //(1:Lattice, 0:Ideal)
    
    int profileType;
    std::string profile_input_file;
    
    double T0;
    
    int fo_type;
    double temp_fo;
    std::string surface_filename_head;
    
    int source;
    


} InitData;


#endif
