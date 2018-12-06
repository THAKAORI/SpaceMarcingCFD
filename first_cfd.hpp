#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlib-cpp-starter/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

class Space_Marching
{
public:
  Space_Marching();
  ~Space_Marching();
  vector<double> boundaryY;
  vector<double> excludeboundary;
  vector<double> momentumboundary;
  vector<double> wallfriction;
  vector<double> coordinatex;
  vector<double> coordinatey;
private:
  double rho;
  double T;
  double viscosity;
  double width;
  double velocity;
  double Re;
  double blasius;
  double nx;
  double dx;
  double ny;
  double ymax;
  double dy;
  double dudy1;
  double dudy2;
  double dudx;
  double dudx1;
  double dudx2;
  double boundaryvelocity;
  double tempboundary;
  double tempboundarymomentum;
  double dpdx;
  vector<int> boundaryindex;
  vector<vector<double> > xvelocitymatrix;
  vector<vector<double> > yvelocitymatrix;
  inline void getviscosity();
  inline void getycoordinate();
  inline void initvelmatrix();
  void update();
  void boundary();
  void wall_friction();
  inline double wall_velocity(int i);
};
