/**
 * @file   step2.cpp
 * @author Wang Heyu <hywang@sixears>
 * @date   Tue Jun  2 17:01:24 2020
 * 
 * @brief  尝试将 AFEPack 对接到我们在 step1 中产生的矩形区域的矩形网格上。
 * 
 * 
 */
#include <iostream>
#include <cmath>
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

#define PI (4.0*atan(1.0))

double u(const double * p)
{
    //return sin(4 * PI * p[0]) * sin(2 * PI * p[1]);
	return sin( 2 * PI * p[0]) * sin(PI * p[1]);
};

double f(const double * p)
{
    //return 20 * PI * PI * u(p);
	return 5 * PI * PI* u( p );
}; 

int main(int argc, char* argv[])
{
    /// 这里基本上和 possion_equation 中配置一致。对比
    /// possion_equation_manual 看更清楚。
    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform;
    rectangle_coord_transform.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof(rectangle_template_geometry);
    /// 一次元。
    rectangle_template_dof.readData("rectangle.2.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function(rectangle_template_dof);
    rectangle_basis_function.readData("rectangle.2.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    template_element.reinit(rectangle_template_geometry,
			    rectangle_template_dof,
			    rectangle_coord_transform,
			    rectangle_basis_function);

    double volume = template_element.volume();
    /// 取了 2 次代数精度。
    const QuadratureInfo<2>& quad_info = template_element.findQuadratureInfo(4);
//    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point = quad_info.quadraturePoint();
//    std::vector<std::vector<std::vector<double> > > basis_gradient = template_element.basisFunction_gradient(q_point);
//    std::vector<std::vector<double> > basis_value = template_element.basisFunction(q_point);
    int n_element_dof = template_element.n_dof();
    int n_bas = rectangle_basis_function.size();

    std::cout << "template element volume:" << volume << std::endl;
    std::cout << "no. of dofs in the element:" << n_element_dof << std::endl;
    std::cout << "no. of quadrature points:" << n_quadrature_point << std::endl;
    /// 产生一个具体单元顶点的缓存。
    double ** arr = (double **) new double* [9];
    for (int i = 0; i < 9; i++)
	arr[i] = (double *) new double [2];
    std::vector<AFEPack::Point<2> > gv(9);
    std::vector<AFEPack::Point<2> > lv(9);

    /// 观察一下模板单元中的自由度、基函数和基函数在具体积分点取值的情
    /// 况。
    for (int i = 0; i < n_element_dof; i++)
    {
	AFEPack::Point<2> pnt = q_point[i];
	/// 第 i 个积分点。
	arr[0][0] = -1.0;arr[0][1] = -1.0;
	arr[1][0] = 1.0;arr[1][1] = -1.0;
	arr[2][0] = 1.0;arr[2][1] = 1.0;
	arr[3][0] = -1.0;arr[3][1] = 1.0; 
	arr[4][0] = 0.0;arr[4][1] = -1.0;
	arr[5][0] = 1.0;arr[5][1] = 0.0;
	arr[6][0] = 0.0;arr[6][1] = 1.0;
	arr[7][0] = -1.0;arr[7][1] = 0.0;
	arr[8][0] = 0.0;arr[8][1] =0.0;
		//for (int j = 0; j < n_bas; j++)
	//{
		
	   // std::cout << "value of basis function " << j << ": " << rectangle_basis_function[j].value(pnt, (const double**)arr) << std::endl;
	//}
    }
    double xa = 0.0;	
    double ya = 0.0;
    double xb = 1.0;
    double yb = 1.0;
    int n = 2;
	int dim = (2 * n + 1) * (2 * n + 1);
	
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
	//for(int i=0;i<dim;i++)
		//std::cout<<boundary[i]<<std::endl;
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

	//FullMatrix<double> stiff_mat(dim,dim);
	//for(int i=0;i<dim;i++)
		//std::cout<<nozeroperow[i]<<std::endl;
	SparsityPattern sp_stiff_matrix(dim, nozeroperow);
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
	    	int idx0 = (2 * j) * (2 * n + 1) + (2 * i) ; 	
	    	int idx1 = (2 * j) * (2 * n + 1) + (2 * i + 2) ;
	    	int idx2 = (2 * j + 2) * (2 * n + 1) + (2 * i + 2) ;
	   		int idx3 = (2 * j + 2) * (2 * n + 1) + (2 * i) ; 
	    	int idx4 = (2 * j) * (2 * n + 1) + (2 * i + 1) ;
			int idx5 = (2 * j + 1) * (2 * n + 1) + (2 * i + 2) ;
			int idx6 = (2 * j + 2) * (2 * n + 1) + (2 * i + 1) ;
			int idx7 = (2 * j + 1) * (2 * n + 1) + (2 * i) ;
			int idx8 = (2 * j + 1) * (2 * n + 1) + (2 * i + 1) ;
			sp_stiff_matrix.add(idx0,idx0);
			sp_stiff_matrix.add(idx0,idx1);
			sp_stiff_matrix.add(idx0,idx2);
			sp_stiff_matrix.add(idx0,idx3);
			sp_stiff_matrix.add(idx0,idx4);
			sp_stiff_matrix.add(idx0,idx5);
			sp_stiff_matrix.add(idx0,idx6);
			sp_stiff_matrix.add(idx0,idx7);
			sp_stiff_matrix.add(idx0,idx8);

			sp_stiff_matrix.add(idx1,idx0);
			sp_stiff_matrix.add(idx1,idx1);
			sp_stiff_matrix.add(idx1,idx2);
			sp_stiff_matrix.add(idx1,idx3);
			sp_stiff_matrix.add(idx1,idx4);
			sp_stiff_matrix.add(idx1,idx5);
			sp_stiff_matrix.add(idx1,idx6);
			sp_stiff_matrix.add(idx1,idx7);
			sp_stiff_matrix.add(idx1,idx8);

			sp_stiff_matrix.add(idx2,idx0);
			sp_stiff_matrix.add(idx2,idx1);
			sp_stiff_matrix.add(idx2,idx2);
			sp_stiff_matrix.add(idx2,idx3);
			sp_stiff_matrix.add(idx2,idx4);
			sp_stiff_matrix.add(idx2,idx5);
			sp_stiff_matrix.add(idx2,idx6);
			sp_stiff_matrix.add(idx2,idx7);
			sp_stiff_matrix.add(idx2,idx8);

			sp_stiff_matrix.add(idx3,idx0);
			sp_stiff_matrix.add(idx3,idx1);
			sp_stiff_matrix.add(idx3,idx2);
			sp_stiff_matrix.add(idx3,idx3);
			sp_stiff_matrix.add(idx3,idx4);
			sp_stiff_matrix.add(idx3,idx5);
			sp_stiff_matrix.add(idx3,idx6);
			sp_stiff_matrix.add(idx3,idx7);
			sp_stiff_matrix.add(idx3,idx8);

			sp_stiff_matrix.add(idx4,idx0);
			sp_stiff_matrix.add(idx4,idx1);
			sp_stiff_matrix.add(idx4,idx2);
			sp_stiff_matrix.add(idx4,idx3);
			sp_stiff_matrix.add(idx4,idx4);
			sp_stiff_matrix.add(idx4,idx5);
			sp_stiff_matrix.add(idx4,idx6);
			sp_stiff_matrix.add(idx4,idx7);
			sp_stiff_matrix.add(idx4,idx8);

			sp_stiff_matrix.add(idx5,idx0);
			sp_stiff_matrix.add(idx5,idx1);
			sp_stiff_matrix.add(idx5,idx2);
			sp_stiff_matrix.add(idx5,idx3);
			sp_stiff_matrix.add(idx5,idx4);
			sp_stiff_matrix.add(idx5,idx5);
			sp_stiff_matrix.add(idx5,idx6);
			sp_stiff_matrix.add(idx5,idx7);
			sp_stiff_matrix.add(idx5,idx8);

			sp_stiff_matrix.add(idx6,idx0);
			sp_stiff_matrix.add(idx6,idx1);
			sp_stiff_matrix.add(idx6,idx2);
			sp_stiff_matrix.add(idx6,idx3);
			sp_stiff_matrix.add(idx6,idx4);
			sp_stiff_matrix.add(idx6,idx5);
			sp_stiff_matrix.add(idx6,idx6);
			sp_stiff_matrix.add(idx6,idx7);
			sp_stiff_matrix.add(idx6,idx8);

			sp_stiff_matrix.add(idx7,idx0);
			sp_stiff_matrix.add(idx7,idx1);
			sp_stiff_matrix.add(idx7,idx2);
			sp_stiff_matrix.add(idx7,idx3);
			sp_stiff_matrix.add(idx7,idx4);
			sp_stiff_matrix.add(idx7,idx5);
			sp_stiff_matrix.add(idx7,idx6);
			sp_stiff_matrix.add(idx7,idx7);
			sp_stiff_matrix.add(idx7,idx8);

			sp_stiff_matrix.add(idx8,idx0);
			sp_stiff_matrix.add(idx8,idx1);
			sp_stiff_matrix.add(idx8,idx2);
			sp_stiff_matrix.add(idx8,idx3);
			sp_stiff_matrix.add(idx8,idx4);
			sp_stiff_matrix.add(idx8,idx5);
			sp_stiff_matrix.add(idx8,idx6);
			sp_stiff_matrix.add(idx8,idx7);
			sp_stiff_matrix.add(idx8,idx8);
		}

	sp_stiff_matrix.compress();
	
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);

    double h = (xb - xa) / n;
    for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
	    	double x0 = ((2 * n - 2 * i) * xa + 2 * i * xb) / (2 * n);
	    	double y0 = ((2 * n - 2 * j) * ya + 2 * j * yb) / (2 * n);
	    	int idx0 = 2 * j * (2 * n + 1) + 2 * i; 

			double x1 = ((2 * n  - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y1 = ((2 * n  - 2 * j) * ya + 2 * j * yb) / (2 * n);
	    	int idx1 = 2 * j * (2 * n + 1) + 2 * i + 2;

	    	double x2 = ((2 * n - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y2 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
	    	int idx2 = (2 * j + 2) * (2 * n + 1) + 2 * i + 2;

	    	double x3 = ((2 * n - 2 * i) * xa + 2 * i * xb) / (2 * n);
	    	double y3 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
	   		int idx3 = (2 * j + 2) * (2 * n + 1) + 2 * i; 

            double x4 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y4 = ((2 * n - 2 * j) * ya + (2 * j ) * yb) / (2 * n);
	    	int idx4 = 2 * j * (2 * n + 1) + 2 * i + 1 ;

			double x5 = ((2 * n - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y5 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx5 = (2 * j + 1) * (2 * n + 1) + 2 * i + 2;

			double x6 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y6 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
			int idx6 = (2 * j + 2) * (2 * n + 1) + 2 * i + 1;

			double x7 = ((2 * n - 2 * i) * xa + (2 * i) * xb) / (2 * n);
	    	double y7 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx7 = (2 * j + 1) * (2 * n + 1) + 2 * i ;

			double x8 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y8 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx8 = (2 * j + 1) * (2 * n + 1) + 2 * i + 1;
		
	    int ele_idx = j * n + i;

	    gv[0][0] = x0;gv[0][1] = y0;  
	    gv[1][0] = x1;gv[1][1] = y1;
	    gv[2][0] = x2;gv[2][1] = y2;
	    gv[3][0] = x3;gv[3][1] = y3;
	    gv[4][0] = x4;gv[4][1] = y4;
	    gv[5][0] = x5;gv[5][1] = y5;
	    gv[6][0] = x6;gv[6][1] = y6;
	    gv[7][0] = x7;gv[7][1] = y7;
	    gv[8][0] = x8;gv[8][1] = y8;

	    lv[0][0] = arr[0][0];lv[0][1] = arr[0][1];
	    lv[1][0] = arr[1][0];lv[1][1] = arr[1][1];
	   	lv[2][0] = arr[2][0];lv[2][1] = arr[2][1];
	    lv[3][0] = arr[3][0];lv[3][1] = arr[3][1];
	    lv[4][0] = arr[4][0];lv[4][1] = arr[4][1];
	    lv[5][0] = arr[5][0];lv[5][1] = arr[5][1];
	   	lv[6][0] = arr[6][0];lv[6][1] = arr[6][1];
	    lv[7][0] = arr[7][0];lv[7][1] = arr[7][1];
		lv[8][0] = arr[8][0];lv[8][1] = arr[8][1];

	    //std::cout << ele_idx << ": " << std::endl;
		//std::cout << idx0<< "->"<< idx1<< "->"<< idx2<< "->"<< idx3<< "->"<< idx4<< "->"<< idx5<< "->"<< idx6<< "->"<< idx7<< "->"<< idx8<<std::endl;
	    std::cout << idx0 << ":(" << x0 << "," << y0 << ") -> "
	    	      << idx1 << ":(" << x1 << "," << y1 << ") -> "
	    	      << idx2 << ":(" << x2 << "," << y2 << ") -> "
	    	      << idx3 << ":(" << x3 << "," << y3 << ") -> "
				  << idx4 << ":(" << x4 << "," << y4 << ") -> "
	    	      << idx5 << ":(" << x5 << "," << y5 << ") -> "
	    	      << idx6 << ":(" << x6 << "," << y6 << ") -> "
				  << idx7 << ":(" << x7 << "," << y7 << ") -> "
	    	      << idx8 << ":(" << x8 << "," << y8 << ")"<<std::endl;
	    // 现在尝试输出具体每个单元的积分点。
	   for (int l = 0; l < n_quadrature_point; l++)
		{
		auto point=rectangle_coord_transform.local_to_global(q_point, lv, gv);
		//std::cout<<quad_info.weight(l)<<std::endl;
		//std::cout<<rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)<<std::endl;
		double Jxy=quad_info.weight(l)*rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)*volume;
		//std::cout<<Jxy<<std::endl;
	    stiff_mat.add(idx0,idx0,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx1,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx2,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx3,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx4,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx5,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx6,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx7,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx0,idx8,Jxy*innerProduct(rectangle_basis_function[0].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx1,idx0,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx1,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx2,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx3,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx4,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx5,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx6,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx7,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx1,idx8,Jxy*innerProduct(rectangle_basis_function[1].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx2,idx0,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx1,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx2,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx3,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx4,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx5,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx6,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx7,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx2,idx8,Jxy*innerProduct(rectangle_basis_function[2].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx3,idx0,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx1,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx2,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx3,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx4,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx5,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx6,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx7,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx3,idx8,Jxy*innerProduct(rectangle_basis_function[3].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx4,idx0,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx1,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx2,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx3,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx4,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx5,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx6,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx7,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx4,idx8,Jxy*innerProduct(rectangle_basis_function[4].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx5,idx0,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx1,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx2,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx3,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx4,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx5,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx6,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx7,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx5,idx8,Jxy*innerProduct(rectangle_basis_function[5].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));

		stiff_mat.add(idx6,idx0,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx1,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx2,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx3,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx4,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx5,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx6,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx7,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx6,idx8,Jxy*innerProduct(rectangle_basis_function[6].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx7,idx0,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx1,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx2,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx3,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx4,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx5,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx6,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx7,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx7,idx8,Jxy*innerProduct(rectangle_basis_function[7].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		stiff_mat.add(idx8,idx0,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[0].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx1,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[1].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx2,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[2].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx3,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[3].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx4,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[4].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx5,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[5].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx6,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[6].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx7,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[7].gradient(point[l],gv)));
		stiff_mat.add(idx8,idx8,Jxy*innerProduct(rectangle_basis_function[8].gradient(point[l],gv),rectangle_basis_function[8].gradient(point[l],gv)));
		
		rhs(idx0)+=Jxy*f(point[l])*rectangle_basis_function[0].value(point[l],gv);
		rhs(idx1)+=Jxy*f(point[l])*rectangle_basis_function[1].value(point[l],gv);
		rhs(idx2)+=Jxy*f(point[l])*rectangle_basis_function[2].value(point[l],gv);
		rhs(idx3)+=Jxy*f(point[l])*rectangle_basis_function[3].value(point[l],gv);
		rhs(idx4)+=Jxy*f(point[l])*rectangle_basis_function[4].value(point[l],gv);
		rhs(idx5)+=Jxy*f(point[l])*rectangle_basis_function[5].value(point[l],gv);
		rhs(idx6)+=Jxy*f(point[l])*rectangle_basis_function[6].value(point[l],gv);
		rhs(idx7)+=Jxy*f(point[l])*rectangle_basis_function[7].value(point[l],gv);
		rhs(idx8)+=Jxy*f(point[l])*rectangle_basis_function[8].value(point[l],gv);
		}/// TO DO: 计算每个积分点上的基函数梯度值，数值积分，拼装局部刚度矩阵，累加至整体刚度矩阵。
	}

	for(int j=0;j<=2*n;j++)
	for(int i=0;i<=2*n;i++)
	{
		int index = j * ( 2 * n + 1 ) + i;
		if(boundary[index]==1)
		{	
			std::cout<<index<<std::endl;
			double x =((2 * n - i )* xa + i * xb)/( 2 * n);
			double y =((2 * n - j )* ya + j * yb)/( 2 * n);
			//std::cout<<"("<<x<<","<<y<<")"<<std::endl;
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
	}
	///　用代数多重网格(AMG)计算线性方程
	AMGSolver solver(stiff_mat);
    /// 这里设置线性求解器的收敛判定为机器 epsilon 乘以矩阵的阶数，也
    /// 就是自由度总数。这个参数基本上是理论可以达到的极限。
	Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);
	std::ofstream f;
	f.open("stiff_matrix.txt");
	f<<"------------"<<std::endl;
	stiff_mat.print_formatted(f,3,false);
	f<<"------------"<<std::endl;
	for(int i =0 ;i<dim;i++)
		f<<rhs[i]<<std::endl;
	f<<"------------"<<std::endl;
	for(int i =0 ;i<dim;i++)
		f<<solution[i]<<std::endl;	
	std::ofstream fs;
	///　输出到output.m用Matlab或Octave运行，得到计算结果
	fs.open("output.m");
	fs<<"x=0:1/"<<(2*n)<<":1;"<<std::endl;
	fs<<"y=0:1/"<<(2*n)<<":1;"<<std::endl;
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

    return 0;
    /// 边界条件。
    /// 矩阵求解。
};