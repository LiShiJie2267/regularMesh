/**
 * @file RectangleDomain.h
 * @author Lishijie (lsj1018845759@outlook.com)
 * @brief 实现矩形区域初始化，并可以实现网格划分.
 * @version 0.1
 * @date 2020-06-23
 * 
 * @copyright Copyright (c) 2020
 * 
 */
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
/**
 * @brief 矩形区域模板类，模板参数为D，表示矩形区域维度；
 * 
 */
TEMPLATE_RECTANGLEDOMAIN
class RectangleDomain
{
    public:    
    AFEPack::Point<D> bottom_left_point;
    AFEPack::Point<D> top_right_point;
    std::vector<int> Divide_vector;
    std::vector<double> interval_length;
    RectangleDomain_Mesh Mesh;
    std::vector<double> boundary;
    std::string mode;
    public:
    /**
     * @brief Construct a new Rectangle Domain object
     * 
     */
    RectangleDomain();
    /**
     * @brief Construct a new Rectangle Domain object
     * 
     * @param p1 矩形区域左下角点
     * @param p2 矩形区域左下角点
     * @param d  std::vector<int> 划分向量
     */
    RectangleDomain(const AFEPack::Point<D> p1,
                    const AFEPack::Point<D> p2,
                    std::vector<int> d);
    /**
     * @brief Construct a new Rectangle Domain object
     * 
     * @param p1 矩形区域左下角点
     * @param p2 矩形区域左下角点
     * @param d  std::vector<int> 划分向量
     * @param length std::vector<double> 向量
     */
    RectangleDomain(const AFEPack::Point<D> p1,
                    const AFEPack::Point<D> p2,
                    std::vector<int> d,
                    std::vector<double> length);
    /**
     * @brief Destroy the Rectangle Domain object
     * 
     */
    ~RectangleDomain(){};

    /**
     * @brief 一维矩形计算区域初始化，可以直接利用RectangleDomain<1> domain
     *        domain.initial_1D_rectangle_domain(x0,x1,nx)进行初始化；
     * @param x0 右上角点x坐标；
     * @param x1 右上角点x坐标；
     * @param nx x方向上划分段数；
     */
    void initial_1D_rectangle_domain(double x0,double x1,int nx);

    /**
     * @brief 二维矩形计算区域初始化，可以直接利用RectangleDomain<2> domain
     *        domain.initial_2D_rectangle_domain(x0,x1，y0,y1,nx,ny)进行初始化；
     * 
     * @param x0 左下角点x坐标；
     * @param x1 右上角点x坐标；
     * @param y0 左下角点y坐标；
     * @param y1 右上角点y坐标；
     * @param nx x方向上划分段数；
     * @param ny y方向上划分段数；
     */
    void initial_2D_rectangle_domain(double x0,double x1,double y0,double y1,int nx,int ny);
    /**
     * @brief 三维矩形计算区域初始化，可以直接利用RectangleDomain<3> domain
     *        domain.initial_3D_rectangle_domain(x0,x1,y0,y1,z1,z2,nx,ny,nz)进行初始化；
     * 
     * @param x0 左下角点x坐标；
     * @param x1 右上角点x坐标；
     * @param y0 左下角点y坐标；
     * @param y1 右上角点y坐标；
     * @param z0 左下角点z坐标；
     * @param z1 右上角点z坐标；
     * @param nx x方向上划分段数；
     * @param ny y方向上划分段数；
     * @param nz z方向上划分段数；
     */
    void initial_3D_rectangle_domain(double x0,double x1,double y0,double y1,double z0,double z1,int nx,int ny,int nz);
    
    /**
     * @brief Get the bottom left point object
     * 
     * @return AFEPack::Point<D> 
     */
    AFEPack::Point<D> get_bottom_left_point();
    
    /**
     * @brief Get the top right point object
     * 
     * @return AFEPack::Point<D> 
     */
    AFEPack::Point<D> get_top_right_point();

    /**
     * @brief Get the Divide vector object
     * 
     * @return std::vector<int> 
     */
    std::vector<int>  get_Divide_vector();

    /**
     * @brief Get the Divide object of rank
     * 
     * @param rank 维度
     * @return int 
     */
    int get_Divide(int rank);

    /**
     * @brief 返回第rank个element信息；
     * 
     */
    loacl_to_global_map element(int rank);
    
    /**
     * @brief 有限元单元个数;
     * 
     * @return int 
     */
    int num_element();

    /**
     * @brief Get the interval length object
     * 
     * @return std::vector<double> 
     */
    std::vector<double> get_interval_length();

    /**
     * @brief Get the interval length object of rank
     * 
     * @param rank 维度
     * @return double 
     */
    double get_interval_length(int rank);

    /**
     * @brief Get the point coord of rank;
     * 
     * @param rank 节点编号
     * @return std::vector<double> 
     */
    std::vector<double> get_point_coord(int rank);

    /**
     * @brief Get the point coord of rank in dimension dim; 
     * 
     * @param rank 节点编号
     * @param dim 维度
     * @return double 
     */
    double get_point_coord(int rank,int dim);

    /**
     * @brief Get the divide mode object
     * 
     * @return std::string 
     */
    std::string get_divide_mode(){return mode;};

    /**
     * @brief Set the divide mode object
     * 
     * @param m 划分模式
     */
    void set_divide_mode(std::string m){mode=m;}

    /**
     * @brief 返回第rank个单元对应的初始迭代器；
     * 
     * @param rank 
     * @return index_iterator 
     */
    index_iterator begin_of_element(int rank){return this->element(rank).begin();};

    /**
     * @brief 返回第rank个单元对应的末尾迭代器；
     * 
     * @param rank 
     * @return index_iterator 
     */
    index_iterator end_of_element(int rank){return this->element(rank).end();};
    /**
     * @brief Generate the mesh with divide mode;
     * 
     */
    void generate_mesh();

    /**
     * @brief Generate the boundary of this rectangle domain;
     * 
     */
    void generate_boundary();
    
    /**
     * @brief 从文件s中读取矩形区域信息，进行初始化；
     * 
     * @param filename 文件名 
     */
    void read_RectangleDomain_Data(std::string filename);

    /**
     * @brief 重载输出运算符；
     * 
     * @tparam DIM 维度；
     * @param out 输出形式；
     * @param domain 待输出对象；
     * @return std::ostream& 
     */
    template <int DIM>
	friend std::ostream& operator<<(std::ostream &out,RectangleDomain<DIM>& domain);

    /**
     * @brief 重载输入运算符；
     * 
     * @tparam DIM 维度
     * @param in 输入格式
     * @param domain 待输入对象；
     * @return std::istream& 
     */
    template <int DIM>
	friend std::istream& operator>>(std::istream &in,RectangleDomain<DIM>& domain);
};

template <int DIM>
std::ostream& operator<<(std::ostream &out,RectangleDomain<DIM>& domain)
{
    out<<"---------The Rectangle Domain-----------"<<std::endl;
    out<<"The dimension of this rectangle domain is : "<<DIM<<std::endl;
    out<<"                                            "<<std::endl;
    out<<"The divide mode of this rectangle domain is : "<<domain.mode<<std::endl; 
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
        for(int i = 0;i < (2 * domain.Divide_vector[0] + 1) * (2 * domain.Divide_vector[1] + 1);i++)
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

