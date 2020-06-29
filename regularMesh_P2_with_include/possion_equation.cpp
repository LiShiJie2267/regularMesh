/**
 * @file possion_equation.cpp
 * @author Lishijie (lsj1018845759@outlook.com)
 * @brief  P2有限元计算possion方程；
 * @version 0.1
 * @date 2020-06-29
 * 
 * @copyright Copyright (c) 2020
 * 
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
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/sparse_mic.h>
#include <lac/sparse_decomposition.h>
#include <lac/full_matrix.h>
#include "../include/RectangleDomain.h"

#define PI (4.0*atan(1.0))
double u(const double * p)
{
    return sin(2 * PI * p[0]) * sin(PI * p[1]);
};

double f(const double * p)
{
    return 5 * PI * PI * u(p);
}; 

typedef std::unordered_map<unsigned int,int> index_map;
int main(int argc, char* argv[])
{
    TemplateGeometry<2> triangle_template_geometry;
    triangle_template_geometry.readData("triangle.tmp_geo");
    CoordTransform<2, 2> triangle_coord_transform;
    triangle_coord_transform.readData("triangle.crd_trs");
    TemplateDOF<2> triangle_template_dof(triangle_template_geometry);
    /// 一次元。
    triangle_template_dof.readData("triangle.2.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> triangle_basis_function(triangle_template_dof);
    triangle_basis_function.readData("triangle.2.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    template_element.reinit(triangle_template_geometry,
			    triangle_template_dof,
			    triangle_coord_transform,
			    triangle_basis_function);

    double volume = template_element.volume();
    const QuadratureInfo<2>& quad_info = template_element.findQuadratureInfo(4);   
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point = quad_info.quadraturePoint();
    int n_element_dof = template_element.n_dof();
    int n_bas = triangle_basis_function.size();
	std::cout<<n_element_dof<<std::endl;
    double ** arr = (double **) new double* [3];
    for (int i = 0; i < 3; i++)
	arr[i] = (double *) new double [2];
    std::vector<AFEPack::Point<2> > gv(6);
    std::vector<AFEPack::Point<2> > lv(6);
	lv[0][0] = 0.0;lv[0][1] = 0.0;
	lv[1][0] = 1.0;lv[1][1] = 0.0;
	lv[2][0] = 0.0;lv[2][1] = 1.0;
	lv[3][0] = 0.5;lv[3][1] = 0.5;
	lv[4][0] = 0.0;lv[4][1] = 0.5;
	lv[5][0] = 0.5;lv[5][1] = 0.0;
    double x0 = 0.0;	
    double y0 = 0.0;
    double x1 = 1.0;
    double y1 = 1.0;
	/// 设置剖分断数和节点总数
    int n = 30;
	int nx = n;
	int ny = n;
	double h = (x1 - x0)/ n;
	RectangleDomain<2> domain;
	/// 利用网格自动生成函数；
	domain.set_divide_mode("P2");
	domain.initial_2D_rectangle_domain(x0,x1,y0,y1,nx,ny);
	domain.generate_mesh();
	domain.generate_boundary();
	int dim = ( 2 * domain.get_Divide( 0 ) + 1 ) * ( 2 * domain.get_Divide( 1 ) + 1 );
	Vector<double> rhs(dim);
	std::vector<unsigned int> nozeroperow(dim);
	for(int i = 0;i <dim ; i++)
		nozeroperow[i]=30;
	SparsityPattern sp_stiff_matrix(dim, nozeroperow);
	for(element_iterator it_b=domain.Mesh.begin();it_b!=domain.Mesh.end();++it_b)
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
			for(index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
				{	
					sp_stiff_matrix.add(it->second,it2->second);
				}
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
			std::vector<AFEPack::Point<2> > point=triangle_coord_transform.local_to_global(q_point, lv, gv);
			double Jxy=quad_info.weight(l)*triangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv);
			for(index_iterator it=it_b2->begin();it!=it_b2->end();it++)
			{
				for(index_iterator it2=it_b2->begin();it2!=it_b2->end();it2++)	
					stiff_mat.add(it->second,it2->second,Jxy*innerProduct(triangle_basis_function[it->first.second].gradient(point[l],gv),triangle_basis_function[it2->first.second].gradient(point[l],gv)));
				rhs(it->second)+=Jxy*f(point[l])*triangle_basis_function[it->first.second].value(point[l],gv);
			}
		}	
	}
	///处理所有边界条件
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

	///　用代数多重网格(AMG)计算线性方程
	AMGSolver solver(stiff_mat);
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
	}
	fs<<"]"<<std::endl;
	fs<<"surf(x,y,U);"<<std::endl;
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
