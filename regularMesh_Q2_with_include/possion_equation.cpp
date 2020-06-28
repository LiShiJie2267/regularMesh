/**
 * @file   possion_equation.cpp
 * @author Lsj <Lsj@Lsj>
 * @date   Tue Jun  5 9:12:34 2020
 * 
 * @brief  利用AFEPack里的函数以及功能，划分结构化网格，使用lagrange九节点二次基函数
 * 		　　，实现Possion方程编值问题计算．
 */
#include <iostream>
#include <cmath>
#include <unordered_map>
#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/SparseMatrixTool.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include "../include/RectangleDomain.h"

#define PI (4.0*atan(1.0))

double u(const double * p)
{
	return sin( 2 * PI * p[0]) * sin(PI * p[1]);
};

double f(const double * p)
{
	return 5 * PI * PI* u( p );
}; 

int main(int argc, char* argv[])
{
    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform;
    rectangle_coord_transform.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof(rectangle_template_geometry);
    /// 一次元。
    rectangle_template_dof.readData("rectangle.lagrange.2.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function(rectangle_template_dof);
    rectangle_basis_function.readData("rectangle.lagrange.2.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    template_element.reinit(rectangle_template_geometry,
			    rectangle_template_dof,
			    rectangle_coord_transform,
			    rectangle_basis_function);
				
    double volume = template_element.volume();
    const QuadratureInfo<2>& quad_info = template_element.findQuadratureInfo(4);
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point = quad_info.quadraturePoint();
	int n_element_dof = template_element.n_dof();
    int n_bas = rectangle_basis_function.size();
    std::vector<AFEPack::Point<2> > gv(9);
    std::vector<AFEPack::Point<2> > lv(9);
	lv[0][0] = -1.0;lv[0][1] = -1.0;
	lv[1][0] = 1.0;lv[1][1] = -1.0;
	lv[2][0] = 1.0;lv[2][1] = 1.0;
	lv[3][0] = -1.0;lv[3][1] = 1.0;
	lv[4][0] = 0.0;lv[4][1] = -1.0;
	lv[5][0] = 1.0;lv[5][1] = 0.0;
	lv[6][0] = 0.0;lv[6][1] = 1.0;
	lv[7][0] = -1.0;lv[7][1] = 0.0;
	lv[8][0] = 0.0;lv[8][1] = 0.0;
    double x0 = 0.0;	
    double y0 = 0.0;
    double x1 = 1.0;
    double y1 = 1.0;
    int n = 20;
	int dim = (2 * n + 1) * (2 * n + 1);
	int nx = n;
	int ny = n;
	RectangleDomain<2> domain;
	domain.set_divide_mode("Q2");
	domain.initial_2D_rectangle_domain(x0,x1,y0,y1,nx,ny);
	domain.generate_mesh();
	domain.generate_boundary();
	Vector<double> rhs(dim);
	std::vector<unsigned int> boundary(dim);
	boundary[0]=boundary[2 * n]=boundary[2 * n * (2 * n + 1)]=boundary[dim - 1] = 1;
	for(int i = 1;i < 2 * n;i++)
	{
		boundary[i] = 1;
		boundary[2 * n * (2 * n + 1) + i] = 1;
	}
	for(int i = 1;i <2 * n;i++)
	{
		boundary[(2 * n + 1) * i] = 1;
		boundary[(2 * n + 1) * i + 2 * n] = 1;
	}
	std::vector<unsigned int> nozeroperow(dim);
	for(int i = 0;i <dim;i++)
		nozeroperow[i]=15;
	nozeroperow[0]=nozeroperow[2 * n]=nozeroperow[2 * n * (2 * n + 1)]=nozeroperow[dim - 1] = 9;
	for(int i = 1;i < 2 * n;i = i + 2)
	{
		nozeroperow[i] = 9;
		nozeroperow[2 * n * (2 * n + 1) +  i] = 9;
	}
	for(int i = 1;i < 2 * n;i = i + 2)
	{
		nozeroperow[(2 * n + 1) * i] = 9;
		nozeroperow[(2 * n + 1) * i + 2 * n] = 9;
	}
	for(int i = 1 ;i < 2 * n;i = i + 2)
		for(int j = 1;j < 2 * n;j = j + 2)
			nozeroperow[(2 * n + 1) * i + j] = 9;
	for(int i = 2 ;i < 2 * n ;i = i + 2)
		for(int j = 2 ;j < 2 * n ;j = j + 2)
			nozeroperow[i *(2 * n + 1) + j] = dim;

	SparsityPattern sp_stiff_matrix(dim, nozeroperow);
	for(element_iterator it_b=domain.Mesh.begin();it_b!=domain.Mesh.end();++it_b)
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
			for(index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
					sp_stiff_matrix.add(it->second,it2->second);
	sp_stiff_matrix.compress();
	SparseMatrix<double> stiff_mat(sp_stiff_matrix);
	for(element_iterator it_b2=domain.Mesh.begin();it_b2!=domain.Mesh.end();++it_b2)
	{       
		for(index_iterator it=it_b2->begin();it!=it_b2->end();++it)
		{
			gv[it->first.second][0]=it->first.first[0];
			gv[it->first.second][1]=it->first.first[1];
		}
		for (int l = 0; l < n_quadrature_point; l++)
		{
			std::vector<AFEPack::Point<2> > point=rectangle_coord_transform.local_to_global(q_point, lv, gv);
			double Jxy=quad_info.weight(l)*rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv);
			for(index_iterator it=it_b2->begin();it!=it_b2->end();it++)
			{
				for(index_iterator it2=it_b2->begin();it2!=it_b2->end();it2++)	
					stiff_mat.add(it->second,it2->second,Jxy*innerProduct(rectangle_basis_function[it->first.second].gradient(point[l],gv),rectangle_basis_function[it2->first.second].gradient(point[l],gv)));
				rhs(it->second)+=Jxy*f(point[l])*rectangle_basis_function[it->first.second].value(point[l],gv);
			}
		}	
	}
	
	for(auto& index:domain.boundary)
	{
		
		double x = domain.get_point_coord(index,0);
		double y = domain.get_point_coord(index,1);
		SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(index);
    	SparseMatrix<double>::iterator row_end = stiff_mat.end(index);
    	double diag = row_iterator->value();
		AFEPack::Point<2> bnd_point;
		bnd_point[0]=x;
		bnd_point[1]=y;
    	double bnd_value = u(bnd_point);
        rhs(index) = diag * bnd_value;
    	for ( ++row_iterator; row_iterator != row_end; ++row_iterator)
        {
            row_iterator->value() = 0.0;
    		int k = row_iterator->column();
            SparseMatrix<double>::iterator col_iterator = stiff_mat.begin(k);   
            SparseMatrix<double>::iterator col_end = stiff_mat.end(k);
    	    for (++col_iterator; col_iterator != col_end; ++col_iterator)
    			if (col_iterator->column() == index)
    			    break;
    		if (col_iterator == col_end)
    		{
    			std::cerr << "Error!" << std::endl;
    			exit(-1);
    		}
    		rhs(k) -= col_iterator->value() * bnd_value;
    		col_iterator->value() = 0.0;	
        } 
	}

	AMGSolver solver(stiff_mat);
    /// 这里设置线性求解器的收敛判定为机器 epsilon 乘以矩阵的阶数，也
    /// 就是自由度总数。这个参数基本上是理论可以达到的极限。
	Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);
	std::ofstream fs;
	///　输出到output.m用Matlab或Octave运行，得到计算结果
	fs.open("output.m");
	fs<<"x="<<x0<<":1/"<<2 * n<<":"<<x1<<";"<<std::endl;
	fs<<"y="<<y0<<":1/"<<2 * n<<":"<<y1<<";"<<std::endl;
	fs<<"[X,Y]=meshgrid(x,y);"<<std::endl;
	fs<<"U=[";
	for(int j=0;j<=2*n;j++)
	{	
		for(int i=0;i<=2*n;i++)
			fs<<solution[j * (2 * n + 1) + i]<<" , ";
		fs<<";"<<std::endl;
		;
	}
	fs<<"]"<<std::endl;
	fs<<"surf(x,y,U);"<<std::endl;
	double h =(x1 - x0)/n;
	double error = 0.0;
    for(unsigned int index = 0; index < dim; index++)
    {
	int x_num = index % (2 * n + 1);
	int y_num = index / (2 * n + 1);
	double x = (x_num * h)/2;	
	double y = (y_num * h)/2;
	AFEPack::Point<2> pnt;
	pnt[0] = x;
	pnt[1] = y;

	double d = (u(pnt) - solution[index]);
	error += d*d;
    }
    error = std::sqrt(error);
    std::cerr << "\nL2 error = " << error << ", tol = " << tol << std::endl;
    return 0;
};
