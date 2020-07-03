/**
 * @file RectangleDomain.templates.h
 * @author Lishijie (lsj1018845759@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2020-06-23
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "RectangleDomain.h"

#define TEMPLATE template<int D>

typedef std::pair<std::vector<double>,unsigned int > pointcoord;
typedef std::map<pointcoord,unsigned int > loacl_to_global_map;
typedef std::vector<loacl_to_global_map>  RectangleDomain_Mesh;

#define generate_mesh1D(nx,x0,x1,i,local_index){\
        double x = (((nx) - (i)) * x0 + (i) * x1) / (nx);\
        std::vector<double> point;\
        point.push_back(x);\
        pointcoord pc(point,local_index);\
        loacl_global_map[pc]=(i);\
}

#define generate_mesh2D(nx,ny,x0,x1,y0,y1,i,j,local_index){\
		double x = ((nx - (i)) * x0 + (i) * x1) / (nx);\
	    double y = ((ny - (j)) * y0 + (j) * y1) / (ny);\
        std::vector<double> point;\
        point.push_back(x);\
        point.push_back(y);\
		pointcoord pc(point,local_index);\
		loacl_global_map[pc]=(j) * (nx + 1) + (i);\
		}

#define generate_mesh3D(nx,ny,nz,x0,x1,y0,y1,z0,z1,i,j,k,local_index){\
		double x = ((nx - (i)) * x0 + (i) * x1) / (nx);\
	    double y = ((ny - (j)) * y0 + (j) * y1) / (ny);\
        double z = ((nz - (k)) * z0 + (k) * z1) / (nz);\
        std::vector<double> point(3);\
        point[0] = x;\
        point[1] = y;\
        point[2] = z;\
		pointcoord pc(point,local_index);\
		loacl_global_map[pc]=(k) * ((nx + 1) * (ny + 1)) + (j) * (nx + 1) + (i);\
		}

#define get_1Dpointcoord(nx,x0,x1,index){\
    double x = (((nx) - (index)) * x0 + (index) * x1) / (nx);\
    point_coord.push_back(x);\
}
#define get_2Dpointcoord(nx,ny,x0,x1,y0,y1,index){\
    int x_num=index % (nx + 1);\
	int y_num=index / (nx + 1);\
	double x = x_num *(x1 - x0)/(nx);\
	double y = y_num *(y1 - y0)/(ny);\
    point_coord.push_back(x);\
    point_coord.push_back(y);\
}

#define get_3Dpointcoord(nx,ny,nz,x0,x1,y0,y1,z0,z1,index){\
    int x_num=((index) % (((nx) + 1) * ((ny) + 1)))%((nx) + 1);\
	int y_num=((index) % (((nx) + 1) * ((ny) + 1)))/((nx) + 1);\
    int z_num=(index) / (((nx) + 1) * ((ny) + 1));\
	double x = x_num *(x1 - x0)/(nx);\
	double y = y_num *(y1 - y0)/(ny);\
    double z = z_num *(z1 - z0)/(nz);\
    point_coord.push_back(x);\
    point_coord.push_back(y);\
    point_coord.push_back(z);\
}
TEMPLATE

RectangleDomain<D>::RectangleDomain()
{
 switch (D)
    {
        case 1:
            mode="L1";
            break;
        case 2:
            mode="Q1";
            break;
        case 3:
            mode="H1";
            break;
        default:
            break;
    }
};

TEMPLATE
RectangleDomain<D>::RectangleDomain(const AFEPack::Point<D> p1,
                                 const AFEPack::Point<D> p2,
                                 std::vector<int> d,
                                 std::vector<double> length)
                                 {
                                     bottom_left_point = p1;
                                     top_right_point = p2;
                                     Divide_vector=d;
                                     interval_length=length;
                                     switch (D)
                                     {
                                     case 1:
                                            mode="L1";
                                            break;
                                     case 2:
                                            mode="Q1";
                                            break;
                                     case 3:
                                            mode="H1";
                                            break;
                                     default:
                                         break;
                                     }
                                 };

TEMPLATE
RectangleDomain<D>::RectangleDomain(const AFEPack::Point<D> p1,
                                 const AFEPack::Point<D> p2,
                                 std::vector<int> d)
                                 {
                                    bottom_left_point = p1;
                                    top_right_point = p2;      
                                    Divide_vector=d;
                                    for(int i = 0; i < D;i++)
                                    {   
                                        interval_length.push_back((p2[i] - p1[i])/Divide_vector[i]);
                                    }
                                    switch (D)
                                     {
                                     case 1:
                                            mode="L1";
                                            break;
                                     case 2:
                                            mode="Q1";
                                            break;
                                     case 3:
                                            mode="H1";
                                            break;
                                     default:
                                         break;
                                     }
                                 };

TEMPLATE        
void RectangleDomain<D>::initial_1D_rectangle_domain(double x0,double x1,int nx)
{
    bottom_left_point[0]=x0;
    top_right_point[0]=x1;
    Divide_vector.push_back(nx);
    interval_length.push_back((x1-x0)/nx);
}

TEMPLATE
void RectangleDomain<D>::initial_2D_rectangle_domain(double x0,double x1,double y0,double y1,int nx,int ny)
{
    bottom_left_point[0]=x0;
    top_right_point[0]=x1;
    bottom_left_point[1]=y0;
    top_right_point[1]=y1;
    Divide_vector.push_back(nx);
    Divide_vector.push_back(ny);
    interval_length.push_back((x1-x0)/nx);
    interval_length.push_back((y1-y0)/ny);
};
TEMPLATE
void RectangleDomain<D>::initial_3D_rectangle_domain(double x0,double x1,double y0,double y1,double z0,double z1,int nx,int ny,int nz)
{
    bottom_left_point[0]=x0;
    top_right_point[0]=x1;
    bottom_left_point[1]=y0;
    top_right_point[1]=y1;
    bottom_left_point[2]=z0;
    top_right_point[2]=z1;
    Divide_vector.push_back(nx);
    Divide_vector.push_back(ny);
    Divide_vector.push_back(nz);
    interval_length.push_back((x1-x0)/nx);
    interval_length.push_back((y1-y0)/ny);
    interval_length.push_back((z1-z0)/nz);
};
TEMPLATE
AFEPack::Point<D> RectangleDomain<D>::get_bottom_left_point()
{
    return bottom_left_point;
};

TEMPLATE
AFEPack::Point<D> RectangleDomain<D>::get_top_right_point()
{
    return top_right_point;
};

TEMPLATE
std::vector<int> RectangleDomain<D>::get_Divide_vector()
{
    return Divide_vector;
};

TEMPLATE
int RectangleDomain<D>::get_Divide(int rank)
{
    assert(rank < D);
    return Divide_vector[rank];
};

TEMPLATE
int RectangleDomain<D>::num_element()
{
    return Mesh.size();
};

TEMPLATE
loacl_to_global_map RectangleDomain<D>::element(int rank)
{
    return Mesh[rank];
};

TEMPLATE
std::vector<double> RectangleDomain<D>::get_interval_length()
{
    return interval_length;
};

TEMPLATE
double RectangleDomain<D>::get_interval_length(int rank)
{   
    assert(rank < D);
    return interval_length[rank];
};

TEMPLATE
void RectangleDomain<D>::generate_mesh()
{
    std::cout<<"-------Generating The Regular Mesh-------"<<std::endl;
    std::cout<<"                   ...                   "<<std::endl;
    switch (D)
    {
    case 1: 
     if(mode=="L1")
        for(int i = 0;i < Divide_vector[0];i++)
        {   
            loacl_to_global_map loacl_global_map;
            generate_mesh1D(Divide_vector[0],bottom_left_point[0],top_right_point[0],i,0);
            generate_mesh1D(Divide_vector[0],bottom_left_point[0],top_right_point[0],i + 1,1);
            Mesh.push_back(loacl_global_map);
        }
    case 2:
    if(mode=="Q1")
        for(int j = 0;j < Divide_vector[1];j++)
        for(int i = 0;i < Divide_vector[0];i++)
        {
		loacl_to_global_map loacl_global_map;
		generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i,j,0);
		generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i + 1,j,1);
		generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i + 1,j + 1,2);
		generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i,j + 1,3);
		Mesh.push_back(loacl_global_map);
        }
    else if(mode=="Q2")
        for (int j = 0; j < Divide_vector[1]; j++)
		for (int i = 0; i < Divide_vector[0]; i++)
		{
		loacl_to_global_map loacl_global_map;
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],2 * i,2 * j,0);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),2 * j,1);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),(2 * j + 2),2);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],2 * i,(2 * j + 2),3);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),2 * j,4);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),(2 * j + 1),5);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j + 2),6);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],2 * i,(2 * j + 1),7);
		generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j + 1),8);
		Mesh.push_back(loacl_global_map);
		}
    else if(mode=="P1")
        {
        for (int j = 0; j < Divide_vector[1]; j++)
		for (int i = 0; i < Divide_vector[0]; i++)
		{
        loacl_to_global_map loacl_global_map;
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i,j,0);
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i+1,j,1);
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i,j+1,2);
        Mesh.push_back(loacl_global_map);       
        }
        for (int j = 0; j < Divide_vector[1]; j++)
		for (int i = 0; i < Divide_vector[0]; i++)
		{
        loacl_to_global_map loacl_global_map;
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(i + 1),(j + 1),0);
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],i,(j + 1),1);
        generate_mesh2D(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(i + 1),j,2);
        Mesh.push_back(loacl_global_map);       
        }
        }
    else if(mode=="P2")
        {
        for (int j = 0; j < Divide_vector[1]; j++)
		for (int i = 0; i < Divide_vector[0]; i++)
		{
        loacl_to_global_map loacl_global_map;
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],2 * i,2 * j,0);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),2 * j,1);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],2 * i,(2 * j + 2),2);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j + 1),3);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i),(2 * j + 1),4);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j),5);
        Mesh.push_back(loacl_global_map);       
        }
        for (int j = 0; j < Divide_vector[1]; j++)
		for (int i = 0; i < Divide_vector[0]; i++)
		{
        loacl_to_global_map loacl_global_map;
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),(2 * j + 2),0);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i),(2 * j + 2),1);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),(2 * j),2);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j + 1),3);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 2),(2 * j + 1),4);
        generate_mesh2D(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                        top_right_point[1],(2 * i + 1),(2 * j + 2),5);
        Mesh.push_back(loacl_global_map);       
        }
        }
    break;
    case 3:
    if(mode=="H1")
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        int ele_idx=k * Divide_vector[2] * Divide_vector[1] + j * Divide_vector[0] + i;
        loacl_to_global_map loacl_global_map;
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k,3);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k + 1,4);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k + 1,5);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k + 1,6);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k + 1,7);
        Mesh.push_back(loacl_global_map);
    }
    else if(mode=="T1")
    {
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        loacl_to_global_map loacl_global_map;/// 0-1-2-5 型单元
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k + 1,3);
        Mesh.push_back(loacl_global_map);
    }
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        loacl_to_global_map loacl_global_map;/// 0-2-3-7 型单元
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k + 1,3);
        Mesh.push_back(loacl_global_map);
    }
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        loacl_to_global_map loacl_global_map;/// 0-4-5-7 型单元
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k + 1,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k + 1,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k + 1,3);
        Mesh.push_back(loacl_global_map);
    }
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        loacl_to_global_map loacl_global_map;/// 2-5-6-7 型单元
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k + 1,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k + 1,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k + 1,3);
        Mesh.push_back(loacl_global_map);
    }
    for(int k = 0;k < Divide_vector[2];k++)
    for(int j = 0;j < Divide_vector[1];j++)
    for(int i = 0;i < Divide_vector[0];i++)
    {
        loacl_to_global_map loacl_global_map;/// 0-5-2-7 型单元
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j,k,0);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j,k + 1,1);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i + 1,j + 1,k,2);
        generate_mesh3D(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],i,j + 1,k + 1,3);
        Mesh.push_back(loacl_global_map);
    }
    }
    default:
        break;
    }
    std::cout<<"                   ...                   "<<std::endl;
    std::cout<<"--------------Generate Done!-------------"<<std::endl;
};

TEMPLATE
std::vector<double> RectangleDomain<D>::get_point_coord(int rank)
{   
    std::vector<double> point_coord;
    switch(D)
    {
    case 1:
        get_1Dpointcoord(Divide_vector[0],bottom_left_point[0],top_right_point[0],rank);
        break;
    case 2:
        if(mode=="Q1"||mode=="P1")
        {
            get_2Dpointcoord(Divide_vector[0],Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                            top_right_point[1],rank);
        }
        else if(mode=="Q2"||mode=="P2")
        {
            get_2Dpointcoord(2 * Divide_vector[0],2 * Divide_vector[1],bottom_left_point[0],top_right_point[0],bottom_left_point[1],
                            top_right_point[1],rank);
        }
        break;
    case 3:
        get_3Dpointcoord(Divide_vector[0],Divide_vector[1],Divide_vector[2],bottom_left_point[0],top_right_point[0],
                        bottom_left_point[1],top_right_point[1],bottom_left_point[2],top_right_point[2],rank);
    default:
        break;
    }
    return point_coord;
};

TEMPLATE
double RectangleDomain<D>::get_point_coord(int rank,int dim)
{
    return get_point_coord(rank)[dim];
};

TEMPLATE
void RectangleDomain<D>::read_RectangleDomain_Data(std::string filename)
{
    std::fstream input;
    input.open(filename);
    input>>(*this);
};

TEMPLATE
void RectangleDomain<D>::generate_boundary()
{
    switch (D)
    {
    case 1:
        boundary.push_back(0);
        boundary.push_back(Divide_vector[0]);
        break;
    case 2:
    if(mode=="Q1"||mode=="P1")
    {
        for(int i = 0; i <= Divide_vector[0];i++)
        {
            boundary.push_back(i);
            boundary.push_back((Divide_vector[0] + 1) * (Divide_vector[1] + 1)- i - 1);
        }
        for(int j = 1;j < Divide_vector[1];j++)
        {
            boundary.push_back((Divide_vector[0] + 1) * j);
            boundary.push_back((Divide_vector[0] + 1) * j + Divide_vector[0]);
        }
    }
    else if(mode=="Q2"||mode=="P2")
    {
        for(int i = 0; i <= 2 * Divide_vector[0];i++)
        {
            boundary.push_back(i);
            boundary.push_back((2 * Divide_vector[0] + 1) * (2 * Divide_vector[1] + 1)- i - 1);
        }
        for(int j = 1;j < 2 * Divide_vector[1];j++)
        {
            boundary.push_back((2 * Divide_vector[0] + 1) * j);
            boundary.push_back((2 * Divide_vector[0] + 1) * j+ 2 * Divide_vector[0]);
        }
    }        
        break;
    case 3:
    if(mode=="H1"||mode=="T1")
    {
        for(int j = 0 ;j <= Divide_vector[1];j++)
        for(int i = 0 ;i <= Divide_vector[0];i++)
        {
            boundary.push_back(j * (Divide_vector[0] + 1) + i);
            boundary.push_back((Divide_vector[0] + 1)*(Divide_vector[1] + 1)*(Divide_vector[2] + 1)-(j * (Divide_vector[0] + 1) + i)-1);
        }
        for(int k = 1;k < Divide_vector[2];k++)
        for(int i = 0;i <= Divide_vector[0];i++)
        {
            boundary.push_back(k * (Divide_vector[0] + 1)*(Divide_vector[1] + 1) + i);
            boundary.push_back((k + 1) * (Divide_vector[0] + 1)*(Divide_vector[1] + 1) -1 - i);
        }
        for(int k = 1;k < Divide_vector[2];k++)
        for(int j = 1;j < Divide_vector[1];j++)
        {
            boundary.push_back(j*(Divide_vector[0]+1)+k*(Divide_vector[0]+1)*(Divide_vector[1]+1));
            boundary.push_back(j*(Divide_vector[0]+1)+k*(Divide_vector[0]+1)*(Divide_vector[1]+1)+Divide_vector[0]);
        }
    }
        break;
    default:
        break;
    }
}

