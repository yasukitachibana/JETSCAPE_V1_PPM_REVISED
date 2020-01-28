#ifndef GRID_H
#define GRID_H

#include <vector>
#include <array>

#include "JetScapeLogger.h"
using namespace Jetscape;

template<class T>
class GridT
{
 private:
    std::vector<T> grid;
      
    int Nx   = 0;
    int Ny   = 0;
    int Neta = 0;
    
    T& get(int ix, int iy, int ieta) {
        return grid[Nx*(Ny*ieta+iy)+ix];
    }
    
 public:
    
    int nX()   const {return(Nx );  }
    int nY()   const {return(Ny );  }
    int nEta() const {return(Neta );}
    int size() const {return Nx*Ny*Neta;}
    int vector_size() const {return grid.size();}

    T& operator()(const int x, const int y, const int eta) {
        assert(0<=x  ); assert(x  <Nx);
        assert(0<=y  ); assert(y  <Ny);
        assert(0<=eta); assert(eta<Neta);
        return get(x, y, eta);
    }

    const T& operator()(int x, int y, int eta) const {
        assert(0<=x  ); assert(x  <Nx);
        assert(0<=y  ); assert(y  <Ny);
        assert(0<=eta); assert(eta<Neta);
        return get(x, y, eta);
    }

    T& operator()(const int i) {
        assert(0<=i  ); assert(i<Nx*Ny*Neta);
        return grid[i];
    }

    const T& operator()(const int i) const {
        assert(0<=i  ); assert(i<Nx*Ny*Neta);
        return grid[i];
    }

    GridT() = default;
    GridT(int Nx_in, int Ny_in, int Neta_in) {
        Nx   = Nx_in  ;
        Ny   = Ny_in  ;
        Neta = Neta_in;
        grid.resize(Nx*Ny*Neta);
    }
    
    void clear() {
        grid.clear();
        grid.shrink_to_fit();
    }
    
    ~GridT(){
        clear();
    }

};




typedef struct {

    double epsilon;
    //double epsilon_prev;

    double rhob;
    double p;
    
    double T;
    double T_prev;
    
    std::array<double, 5> U;
    
    std::array<std::array<double, 5>, 3> U_surf;
    
    //std::array<double, 5> U_prev;
    
    std::array<double, 4> u;
    //std::array<double, 4> u_prev;
    
} FluidCell;

typedef GridT<FluidCell> SCGrid;

#endif
