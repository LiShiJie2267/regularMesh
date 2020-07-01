/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Seven
 * @Date: 2020-06-14 21:07:15
 * @LastEditors: Seven
 * @LastEditTime: 2020-07-01 19:54:38
 */ 

#include "RectangleDomain.h"
#include <vector>
#include <iostream>
int main()
{
    RectangleDomain<1> mesh1D;
    mesh1D.initial_1D_rectangle_domain(0.0,2.0,10);
    mesh1D.generate_mesh();
    std::cout<<mesh1D<<std::endl;

    RectangleDomain<3> mesh3D;
    mesh3D.initial_3D_rectangle_domain(0.0,1.0,0.0,1.0,0.0,1.0,3,3,2);
    //mesh3D.set_divide_mode("P2");
    mesh3D.generate_mesh();
    mesh3D.generate_boundary();
    std::cout<<mesh3D<<std::endl;
}