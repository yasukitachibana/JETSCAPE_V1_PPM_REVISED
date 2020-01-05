#ifndef PATON_ADS_H_
#define PATON_ADS_H_

#include <string>
#include <vector>

class PartonAdS{
private:
    void LoadPartonList(std::string source_input_file);
    void Message();
    double totalE, totalPx, totalPy, totalPz;
    
public:
    PartonAdS(std::string source_input_file);
    ~PartonAdS();// destructor
    

    struct Droplet{
        int sn, pid, status;
        double tau, x, y, eta_s;
        double e, px, py, pz;
        double de, dpx, dpy, dpz;
    };
    std::vector <Droplet> droplets;
    
};

#endif 
