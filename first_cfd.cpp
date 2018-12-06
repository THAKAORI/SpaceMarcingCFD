#include "first_cfd.hpp"

using namespace std;
Space_Marching::Space_Marching()
{
  rho = 1.225;
  T = 293.0;
  Space_Marching::getviscosity();
  width = 0.5;
  velocity = 20.0;
  boundaryvelocity = velocity * 0.995;
  Re = rho * velocity * width / viscosity;
  blasius = width * 5.3 / sqrt(Re);
  nx = 50000;
  dx = width / (double)nx;
  //printf("%f\n",dx);
  ny = 200;
  tempboundary = 0;
  ymax = blasius * 2.0;
  dpdx = 0;
  dy = ymax / (double)ny;
  coordinatex = vector<double>(nx, 0);
  coordinatey = vector<double>(ny, 0);
  boundaryindex = vector<int>(nx, 0);
  boundaryY = vector<double>(nx, 0);
  excludeboundary = vector<double>(nx, 0);
  momentumboundary = vector<double>(nx, 0);
  wallfriction = vector<double>(nx, 0);
  Space_Marching::getycoordinate();
  xvelocitymatrix = vector<vector<double> >(ny, vector<double>(nx + 1, 0));
  yvelocitymatrix = vector<vector<double> >(ny, vector<double>(nx + 1, 0));
  Space_Marching::initvelmatrix();
  printf("initialize ok\n");
  Space_Marching::update();
  printf("update ok\n");
  Space_Marching::boundary();
  printf("boundary ok\n");
  Space_Marching::wall_friction();

  //cout << dy << endl;
}
Space_Marching::~Space_Marching()
{
}
inline void Space_Marching::getviscosity()
{
  viscosity = 1.458 / 1000000 * pow(T, 1.5) / (T + 110.4);
}
inline void Space_Marching::getycoordinate()
{
  for (int i = 0; i < ny; i++)
  {
    coordinatey[i] = dy * i;
  }
  for (int i = 0; i < nx; i++)
  {
    coordinatex[i] = dx * i;
  }
}
inline void Space_Marching::initvelmatrix()
{
  for (int i = 0; i < ny; i++)
  {
    xvelocitymatrix[i][0] = velocity;
    yvelocitymatrix[i][0] = 0.0;
  }
}
void Space_Marching::update()
{
  for (int j = 0; j < nx; j++)
  {
    for (int i = 1; i < ny - 1; i++)
    {
      dudy1 = (xvelocitymatrix[i + 1][j] - xvelocitymatrix[i - 1][j]) / (2 * dy);
      dudy2 = (xvelocitymatrix[i + 1][j] - 2 * xvelocitymatrix[i][j] + xvelocitymatrix[i - 1][j]) / pow(dy, 2);
      dudx = (dudy2 * viscosity / rho - dpdx / rho - yvelocitymatrix[i][j] * dudy1) / xvelocitymatrix[i][j];
      xvelocitymatrix[i][j + 1] = xvelocitymatrix[i][j] + dx * dudx;
      xvelocitymatrix[0][j + 1] = 0;
      xvelocitymatrix[ny - 1][j + 1] = velocity;//xvelocitymatrix[ny - 2][j + 1];
      yvelocitymatrix[0][j + 1] = wall_velocity(j + 1);
    }

    for (int i = 0; i < ny; i++)
    {
      xvelocitymatrix[i][nx] = xvelocitymatrix[i][nx - 1];
    }
    for (int i = 0; i < ny - 1; i++)
    {
      dudx1 = (xvelocitymatrix[i][j + 1] - xvelocitymatrix[i][j]) / dx;
      dudx2 = (xvelocitymatrix[i + 1][j + 1] - xvelocitymatrix[i + 1][j]) / dx;
      yvelocitymatrix[i + 1][j + 1] = yvelocitymatrix[i][j + 1] - (dudx1 + dudx2) * dy / 2.0;
    }
  }
}
inline double Space_Marching::wall_velocity(int i)
{
  return 0;
}
void Space_Marching::boundary()
{
  for (int j = 0; j < nx; j++)
  {
    for (int i = 0; i < ny; i++)
    {
      if (xvelocitymatrix[i][j] > boundaryvelocity)
      {
        boundaryindex[j] = i;
        break;
      }
    }
  }

  for (int j = 0; j < nx; j++)
  {
    boundaryY[j] = coordinatey[boundaryindex[j]];
  }
  printf("nomal boundary ok\n");
  for (int j = 0; j < nx; j++)
  {
    tempboundary = 0;
    tempboundarymomentum = 0;
    for (int i = 0; i < ny; i++)
    {
      tempboundary += (1 - xvelocitymatrix[i][j] / xvelocitymatrix[ny - 1][j]) * dy;
      tempboundarymomentum += xvelocitymatrix[i][j] / xvelocitymatrix[ny - 1][j] * (1 - xvelocitymatrix[i][j] / xvelocitymatrix[ny - 1][j]) * dy;
    }
    //printf("%f\n",xvelocitymatrix[ny - 2][j] );
    excludeboundary[j] = tempboundary;
    momentumboundary[j] = tempboundarymomentum;
  }
  printf("other-boundary ok\n");
  /*plt::title("boundary");
  plt::xlabel("x");
  plt::ylabel("y");
  plt::named_plot("nom_boundary", coordinatex, boundaryY,"-r");
  plt::named_plot("ex_boundary", coordinatex, excludeboundary,"-b");
  plt::named_plot("mom_boundary", coordinatex, momentumboundary,"-y");
  plt::legend();
  plt::save("./boundary.pdf");
  */
}
void Space_Marching::wall_friction()
{
  for (int i = 0; i < nx; i++)
  {
    wallfriction[i] = viscosity * (xvelocitymatrix[1][i] - xvelocitymatrix[0][i]) / dy * 2 / rho / pow(xvelocitymatrix[ny - 1][i], 2);
  }
  /*
  plt::title("friction-coeficient");
  plt::xlabel("x");
  plt::ylabel("friction-coeficient");
  plt::named_plot("tau_w", coordinatex, wallfriction,"-c");
  plt::legend();
  plt::save("./wallfriction.pdf");
  */
}

void graphplotbl(Space_Marching cfd)
{
  plt::title("boundary");
  plt::xlabel("x");
  plt::ylabel("y");
  plt::named_plot("nom_boundary", cfd.coordinatex, cfd.boundaryY,"-r");
  plt::named_plot("ex_boundary", cfd.coordinatex, cfd.excludeboundary,"-b");
  plt::named_plot("mom_boundary", cfd.coordinatex, cfd.momentumboundary,"-y");
  plt::legend();
  plt::save("./boundary.pdf");
}
void graphplotfr(Space_Marching cfd)
{
  plt::title("friction-coeficient");
  plt::xlabel("x");
  plt::ylabel("friction-coeficient");
  plt::named_plot("tau_w", cfd.coordinatex, cfd.wallfriction,"-c");
  plt::legend();
  plt::save("./wallfriction.pdf");
}
int main(int argc, char** argv)
{
  Space_Marching cfd;
  //graphplotbl(cfd);
  graphplotfr(cfd);
  return 0;
}
