#include <math.h>

#include "PPMpartonCloud.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace std;

PartonCloud::PartonCloud(string source_input_file){
    
    JSINFO << "<-[PPM] Creating PartonCloud ->";
    
    totalE = 0.0;
    totalPx=0.0;
    totalPy=0.0;
    totalPz=0.0;
    
    LoadPartonList(source_input_file);
    
}

PartonCloud::~PartonCloud(){
    Message();
    drops.clear();
    JSINFO << "<-[PPM] Deleting PartonCloud ->";
}// destructor

void PartonCloud::LoadPartonList(string source_input_file){

    ifstream source_ifs( source_input_file.c_str() );


    JSINFO
    << "<-[PPM] PartonCloud, Loading Source Term Input File: "
    << source_input_file
    << " ->";
    
    if(source_ifs.fail()) {
        JSINFO << "<-[PPM] Source Term Input File - Found ->";
    }else{
        JSINFO<< "<-[PPM] Source Term Input File - Not Found ->";
    }

    int id_d, pstat_d, dstat_d;
    double  t_d, x_d, y_d, z_d;
    double  e_d, px_d, py_d, pz_d;

    string source_input_line;
    while(getline(source_ifs, source_input_line)){

//        JSINFO
//        << "[ORIGINAL] "
//        << source_input_line;

        dstat_d = 1;
        sscanf(source_input_line.data(),
               "%d %d %lf %lf %lf %lf %lf %lf %lf %lf",
               &id_d, &pstat_d,
               &t_d, &x_d, &y_d, &z_d,
               &e_d, &px_d, &py_d, &pz_d);

        Droplet this_drop;

        this_drop.index = id_d;
        this_drop.pstat = pstat_d;
        this_drop.dstat = dstat_d;
        this_drop.x[0] = t_d;
        this_drop.x[1] = x_d;
        this_drop.x[2] = y_d;
        this_drop.x[3] = z_d;
        this_drop.p[0] = e_d;
        this_drop.p[1] = px_d;
        this_drop.p[2] = py_d;
        this_drop.p[3] = pz_d;

        drops.push_back(this_drop);

//        JSINFO
//        << "[ COPIED ] "
//        << drops[id_d].index << " "
//        << drops[id_d].pstat << " "
//        << drops[id_d].x[0] << " "
//        << drops[id_d].x[1] << " "
//        << drops[id_d].x[2] << " "
//        << drops[id_d].x[3] << " "
//        << drops[id_d].p[0] << " "
//        << drops[id_d].p[1] << " "
//        << drops[id_d].p[2] << " "
//        << drops[id_d].p[3] << " "
//        << drops[id_d].dstat;
        
        if(pstat_d!=-1){
            totalE += e_d;
            totalPx += px_d;
            totalPy += py_d;
            totalPz += pz_d;
        }else{
            totalE -= e_d;
            totalPx -= px_d;
            totalPy -= py_d;
            totalPz -= pz_d;
        }
    }

    source_ifs.close();
    
    Message();

}

void PartonCloud::Message(){
    
    JSINFO
    << "<-[PPM] Source Term Loaded, Total Source Number = "
    << drops.size()
    << " ->";
    
    JSINFO
    << "<-[PPM] Source to be deposited into the fluid: ->";
    JSINFO
    << "<-[PPM] total_e: " << totalE
    << " GeV, total_px: " << totalPx
    << " GeV, total_py: " << totalPy
    << " GeV, total_pz: " << totalPz
    << " GeV, total_mass^2: " << (totalE*totalE - totalPx*totalPx - totalPy*totalPy - totalPz*totalPz)
    << " GeV^2 ->";
    
}
