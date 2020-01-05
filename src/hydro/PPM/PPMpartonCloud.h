#ifndef PATON_CLOUD_H_
#define PATON_CLOUD_H_

#include <string>
#include <vector>

#include "PPMutil.h"


class PartonCloud{
private:
    void LoadPartonList(std::string source_input_file);
    void Message();
    double totalE, totalPx, totalPy, totalPz;
    
public:
    PartonCloud(std::string source_input_file);
    ~PartonCloud();// destructor
    
    public:
    struct Droplet{
        int index;
        int pstat;
        int dstat;
        double x[4];//x,y,z
        double p[4];//E,px,py,pz
    };
    std::vector <Droplet> drops;
    
};

#endif 
