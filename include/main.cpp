/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Seven
 * @Date: 2020-06-14 21:07:15
 * @LastEditors: Seven
 * @LastEditTime: 2020-06-29 19:40:43
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

    RectangleDomain<2> mesh2D;
    mesh2D.initial_2D_rectangle_domain(0.0,1.0,0.0,1.0,2,2);
    mesh2D.set_divide_mode("P2");
    mesh2D.generate_mesh();
    mesh2D.generate_boundary();
    std::cout<<mesh2D<<std::endl;
}