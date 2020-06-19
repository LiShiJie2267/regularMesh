
#include "RectangleDomain.h"
#include <vector>
#include <iostream>
int main()
{
    RectangleDomain<1> mesh1D;
    mesh1D.initial_1D_rectangle_domain(0.0,2.0,10);
    mesh1D.generate_mesh();
    std::cout<<mesh1D<<std::endl;

    RectangleDomain<2> mesh2D;
    mesh2D.initial_2D_rectangle_domain(0.0,1.0,0.0,1.0,2,4);
    mesh2D.generate_mesh();
    std::cout<<mesh2D<<std::endl;

    RectangleDomain<3> mesh3D;
    mesh3D.initial_3D_rectangle_domain(0.0,1.0,0.0,1.0,0.0,1.0,2,2,2);
    mesh3D.generate_mesh();
    std::cout<<mesh3D<<std::endl;

    RectangleDomain<2> m;
    m.read_RectangleDomain_Data("input.txt");
    std::cout<<m;
}