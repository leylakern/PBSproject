#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>
#include <random>

class sph {
public:
    sph();
  //  int pppp;
    Eigen::MatrixX3d x;
    Eigen::MatrixX3d u;
    Eigen::MatrixXi neighbours;
    Eigen::VectorXd rho;
    Eigen::VectorXd p;
    Eigen::MatrixX3d f;
    
    Eigen::MatrixX3d c; //cells
    Eigen::VectorXd lscl; //hold index of atom j to which atom i points
    Eigen::VectorXd head; //hold index of first atom in the c-th cell or (-1) if empty cell

    void render();
    void update();
    
    void step();
    void searchNeighbour(int i);
    double w(double x_ij);
    double w_grad(double x_ij);
    double w_grad2(double x_ij);
    void computeDensity(int i);
    void computePressure(int i);
    void computeForce(int i);
    
    
};
