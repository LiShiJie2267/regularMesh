#ifndef _RECTANGLEDOMAIN_H_
#define _RECTANGLEDOMAIN_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <AFEPack/Geometry.h>


#define TEMPLATE_RECTANGLEDOMAIN template<int D>

typedef std::pair<std::vector<double>,unsigned int > pointcoord;
typedef std::map<pointcoord,unsigned int > loacl_to_global_map;
typedef loacl_to_global_map::iterator index_iterator;
typedef std::vector<loacl_to_global_map>  RectangleDomain_Mesh;
typedef std::vector<loacl_to_global_map>::iterator element_iterator;

TEMPLATE_RECTANGLEDOMAIN
class RectangleDomain
{
    public:    
    AFEPack::Point<D> bottom_left_point;
    AFEPack::Point<D> top_right_point;
    std::vector<int> Divide_vector;
    std::vector<double> interval_length;
    RectangleDomain_Mesh Mesh;
    std::string mode;
    public:

    RectangleDomain();
    RectangleDomain(const AFEPack::Point<D> p1,
                    const AFEPack::Point<D> p2,
                    std::vector<int> d);
                    
    RectangleDomain(const AFEPack::Point<D> p1,
                    const AFEPack::Point<D> p2,
                    std::vector<int> d,
                    std::vector<double> length);

    ~RectangleDomain(){};
    void initial_1D_rectangle_domain(double x0,double x1,int nx);
    void initial_2D_rectangle_domain(double x0,double x1,double y0,double y1,int nx,int ny);
    void initial_3D_rectangle_domain(double x0,double x1,double y0,double y1,double z0,double z1,int nx,int ny,int nz);
    
    AFEPack::Point<D> get_bottom_left_point();
    AFEPack::Point<D> get_top_right_point();
    std::vector<int>  get_Divide_vector();
    int get_Divide(int rank);
    std::vector<double> get_interval_length();
    double get_interval_length(int rank);
    std::vector<double> get_point_coord(int rank);
    double get_point_coord(int rank,int dim);
    std::string get_divide_mode(){return mode;};
    void set_divide_mode(std::string m){mode=m;}
    void generate_mesh();
    void read_RectangleDomain_Data(std::string s);

    template <int DIM>
	friend std::ostream& operator<<(std::ostream &out,RectangleDomain<DIM>& domain);

    template <int DIM>
	friend std::istream& operator>>(std::istream &in,RectangleDomain<DIM>& domain);
};

template <int DIM>
std::ostream& operator<<(std::ostream &out,RectangleDomain<DIM>& domain)
{
    out<<"---------The Rectangle Domain-----------"<<std::endl;
    out<<"The dimension of this rectangle domain is : "<<DIM<<std::endl;
    out<<"                                            "<<std::endl;
    switch (DIM)
    {
    case 1:
        out<<"The interval of x direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[0]<<" , "<< domain.top_right_point[0]<<" ] "<<std::endl;
        out<<"The divison number of x direction is : "<< domain.Divide_vector[0]<<std::endl;
        out<<"The length of each interval in x direction is : "<< domain.interval_length[0]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"---The point coord in this rectangle domain mesh---"<<std::endl;
        for(int i = 0;i < domain.Divide_vector[0] + 1;i++)
        {
            std::cout<<"This point coord of "<<i<<" th is :";
            std::cout<<" ( "<<domain.get_point_coord(i,0)<<" ) "<<std::endl;
        }    
        break;
    case 2:
        out<<"The interval of x direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[0]<<" , "<< domain.top_right_point[0]<<" ] "<<std::endl;
        out<<"The divison number of x direction is : "<< domain.Divide_vector[0]<<std::endl;
        out<<"The length of each interval in x direction is : "<< domain.interval_length[0]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"The interval of y direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[1]<<" , "<< domain.top_right_point[1]<<" ] "<<std::endl;
        out<<"The divison number of y direction is : "<< domain.Divide_vector[1]<<std::endl;
        out<<"The length of each interval in y direction is : "<< domain.interval_length[1]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"---The point coord in this rectangle domain mesh---"<<std::endl;
        for(int i = 0;i < (domain.Divide_vector[0] + 1) * (domain.Divide_vector[1] + 1);i++)
        {   
            std::cout<<"This point coord of "<<i<<" th is :";
            std::cout<<" ( "<<domain.get_point_coord(i,0)<<" , "<<domain.get_point_coord(i,1)<<" ) "<<std::endl;
        }
        break;
    case 3:
        out<<"The interval of x direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[0]<<" , "<< domain.top_right_point[0]<<" ] "<<std::endl;
        out<<"The divison number of x direction is : "<< domain.Divide_vector[0]<<std::endl;
        out<<"The length of each interval in x direction is : "<< domain.interval_length[0]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"The interval of y direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[1]<<" , "<< domain.top_right_point[1]<<" ] "<<std::endl;
        out<<"The divison number of y direction is : "<< domain.Divide_vector[1]<<std::endl;
        out<<"The length of each interval in y direction is : "<< domain.interval_length[1]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"The interval of z direction in this rectangle domain is :";
        out<<" [ "<<domain.bottom_left_point[2]<<" , "<< domain.top_right_point[2]<<" ] "<<std::endl;
        out<<"The divison number of z direction is : "<< domain.Divide_vector[2]<<std::endl;
        out<<"The length of each interval in z direction is : "<< domain.interval_length[2]<<std::endl;
        out<<"                                              "<<std::endl;
        out<<"---The point coord in this rectangle domain mesh---"<<std::endl;
        for(int i = 0;i < (domain.Divide_vector[0] + 1) * (domain.Divide_vector[1] + 1) *(domain.Divide_vector[2] + 1);i++)
        {   
            std::cout<<"This point coord of "<<i<<"th is :";
            std::cout<<" ( "<<domain.get_point_coord(i,0)<<" , "<<domain.get_point_coord(i,1)<<" , "<<domain.get_point_coord(i,2)<<" ) "<<std::endl;
        }
        break;
    default:
        break;
    }
    return out;
};
template <int DIM>
std::istream& operator>>(std::istream &in,RectangleDomain<DIM>& domain)
{
    std::cout<<"---------The Rectangle Domain input---------"<<std::endl;
    std::cout<<"The dimension of this rectangle domain is : "<<DIM<<std::endl;
    //std::cout<<"                                            "<<std::endl;
    //int n;
    //in>>n;
    //DIM=n;
    switch (DIM)
    {
    case 1:
        std::cout<<"The interval of x direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[0];
        in>>domain.top_right_point[0];
        std::cout<<"The divison number of x direction is : "<<std::endl;
        int nx;
        in>>nx;
        domain.Divide_vector.push_back(nx);
        std::cout<<"The length of each interval in x direction is : "<<std::endl;
        double lx;
        in>>lx;
        domain.interval_length.push_back(lx);
        break;
    case 2:
        std::cout<<"The interval of x direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[0];
        in>>domain.top_right_point[0];
        std::cout<<"The divison number of x direction is : "<<std::endl;
        int nx2;
        in>>nx2;
        domain.Divide_vector.push_back(nx2);
        std::cout<<"The length of each interval in x direction is : "<<std::endl;
        double lx2;
        in>>lx2;
        domain.interval_length.push_back(lx2);
        std::cout<<"                                            "<<std::endl;
        std::cout<<"The interval of y direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[1];
        in>>domain.top_right_point[1];
        std::cout<<"The divison number of y direction is : "<<std::endl;
        int ny2;
        in>>ny2;
        domain.Divide_vector.push_back(ny2);
        std::cout<<"The length of each interval in y direction is : "<<std::endl;
        double ly2;
        in>>ly2;
        domain.interval_length.push_back(ly2);
        break;
    case 3:
        std::cout<<"The interval of x direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[0];
        in>>domain.top_right_point[0];
        std::cout<<"The divison number of x direction is : "<<std::endl;
        int nx3;
        in>>nx3;
        domain.Divide_vector.push_back(nx3);
        std::cout<<"The length of each interval in x direction is : "<<std::endl;
        double lx3;
        in>>lx3;
        domain.interval_length.push_back(lx3);

        std::cout<<"                                            "<<std::endl;
        std::cout<<"The interval of y direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[1];
        in>>domain.top_right_point[1];
        std::cout<<"The divison number of y direction is : "<<std::endl;
        int ny3;
        in>>ny3;
        domain.Divide_vector.push_back(ny3);
        std::cout<<"The length of each interval in y direction is : "<<std::endl;
        double ly3;
        in>>ly3;
        domain.interval_length.push_back(ly3);
        std::cout<<"                                            "<<std::endl;

        std::cout<<"The interval of z direction in this rectangle domain is :"<<std::endl;
        in>>domain.bottom_left_point[2];
        in>>domain.top_right_point[2];
        std::cout<<"The divison number of z direction is : "<<std::endl;
        int nz3;
        in>>nz3;
        domain.Divide_vector.push_back(nz3);
        std::cout<<"The length of each interval in z direction is : "<<std::endl;
        double lz3;
        in>>lz3;
        domain.interval_length.push_back(lz3);
        break;
    default:
        break;
    }
    return in;
};
#else
#endif

