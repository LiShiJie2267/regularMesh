/**
 * @file   possion_equation.cpp
 * @author Li ShiJie <lsj@lsj>
 * @date   Tue Jun  10 12:34:11 2020
 * 
 * @brief  利用AFEPack里的函数以及功能，划分结构化网格，使用三角有限元三节点线性基函数
 * 		　　，实现Possion方程编值问题计算．
 * 
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
    triangle_template_dof.readData("triangle.1.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> triangle_basis_function(triangle_template_dof);
    triangle_basis_function.readData("triangle.1.bas_fun");
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

    double ** arr = (double **) new double* [3];
    for (int i = 0; i < 3; i++)
	arr[i] = (double *) new double [2];
    std::vector<AFEPack::Point<2> > gv(3);
	std::vector<AFEPack::Point<2> > gv_2(3);
    std::vector<AFEPack::Point<2> > lv(3);

	arr[0][0] = 0.0;
	arr[0][1] = 0.0;
	arr[1][0] = 1.0;
	arr[1][1] = 0.0;
	arr[2][0] = 0.0;
	arr[2][1] = 1.0;

    double x0 = 0.0;	
    double y0 = 0.0;
    double x1 = 1.0;
    double y1 = 1.0;
	/// 设置剖分断数和节点总数
    int n = 50;
	int dim=(n+1)*(n+1);

	Vector<double> rhs(dim);
	std::vector<unsigned int> nozeroperow(dim);
	for(int i = 0;i <dim ; i++)
		nozeroperow[i]=5;
	nozeroperow[0]=3;
	nozeroperow[dim-1]=3;
	nozeroperow[n]=3;
	nozeroperow[dim-1-n]=3;
	for(int i=1;i<n;i++)
	{
		nozeroperow[i]=4;
		nozeroperow[dim-1-i]=4;
		nozeroperow[i*(n+1)]=4;
		nozeroperow[(i+1)*n+i]=4;
	}
	SparsityPattern sp_stiff_matrix(dim, nozeroperow);

	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
	    int idx00 = j * (n + 1) + i; 
	    int idx10 = j * (n + 1) + i + 1; 
	    int idx01 = (j + 1) * (n + 1) + i; 

		index_map index({{0,idx00},{1,idx10},{2,idx01}},template_element.n_dof());
		for(int k1 = 0;k1 < template_element.n_dof(); k1++)
		{
			for(int k2 = 0;k2 + k1 < template_element.n_dof(); k2++)
				sp_stiff_matrix.add(index[k1],index[k2]);
		}
	}
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
	    int idx10 = j * (n + 1) + i + 1; 
	    int idx11 = (j + 1) * (n + 1) + i + 1; 
	    int idx01 = (j + 1) * (n + 1) + i; 

		index_map index2({{0,idx11},{1,idx01},{2,idx10}},template_element.n_dof());
		for(int k1 = 0;k1 < template_element.n_dof(); k1++)
		{
			for(int k2 = 0;k2 + k1 < template_element.n_dof(); k2++)
				sp_stiff_matrix.add(index2[k1],index2[k2]);
		}
	}
	
	sp_stiff_matrix.compress();
	SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    double h = (x1 - x0) / n;
    for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
	    double x00 = ((n - i) * x0 + i * x1) / n;
	    double y00 = ((n - j) * y0 + j * y1) / n;
	    int idx00 = j * (n + 1) + i; 
	    double x10 = ((n - i - 1) * x0 + (i + 1) * x1) / n;
	    double y10 = ((n - j ) * y0 + j * y1) / n;
	    int idx10 = j * (n + 1) + i + 1; 
	    double x11 = ((n - i - 1) * x0 + (i + 1) * x1) / n;
	    double y11 = ((n - j - 1) * y0 + (j + 1) * y1) / n;
	    int idx11 = (j + 1) * (n + 1) + i + 1; 
	    double x01 = ((n - i) * x0 + i * x1) / n;
	    double y01 = ((n - j - 1) * y0 + (j + 1) * y1) / n;
	    int idx01 = (j + 1) * (n + 1) + i; 
	    
	    int ele_idx = j * n + i;

	    gv[0][0] = x00;
	    gv[0][1] = y00;
	    gv[1][0] = x10;
	    gv[1][1] = y10;
	    gv[2][0] = x01;
	    gv[2][1] = y01;
		gv_2[0][0] = x11;
		gv_2[0][1] = y11;
		gv_2[1][0] = x01;
	    gv_2[1][1] = y01;
	    gv_2[2][0] = x10;
	    gv_2[2][1] = y10;
	    lv[0][0] = arr[0][0];
	    lv[0][1] = arr[0][1];
	    lv[1][0] = arr[1][0];
	    lv[1][1] = arr[1][1];
	    lv[2][0] = arr[2][0];
	    lv[2][1] = arr[2][1];
	    for (int l = 0; l < n_quadrature_point; l++)
		{
		auto point=triangle_coord_transform.local_to_global(q_point, lv, gv);
		auto point_2=triangle_coord_transform.local_to_global(q_point, lv, gv_2);
		double Jxy=quad_info.weight(l)*triangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)*volume;
		double Jxy_2=quad_info.weight(l)*triangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv_2)*volume;
		index_map index3({{0,idx00},{1,idx10},{2,idx01}},template_element.n_dof());
		for(int base1 = 0;base1 < template_element.n_dof(); base1++)
		{
			for(int base2 = 0;base2 < template_element.n_dof(); base2++)
				stiff_mat.add(index3[base1],index3[base2],Jxy*innerProduct(triangle_basis_function[base1].gradient(point[l],gv),triangle_basis_function[base2].gradient(point[l],gv)));
			rhs(index3[base1])+=Jxy*f(point[l])*triangle_basis_function[base1].value(point[l],gv);
		}
		index_map index4({{0,idx11},{1,idx01},{2,idx10}},template_element.n_dof());
		for(int base1 = 0;base1 < template_element.n_dof(); base1++)
		{
			for(int base2 = 0;base2 < template_element.n_dof(); base2++)
				stiff_mat.add(index4[base1],index4[base2],Jxy_2*innerProduct(triangle_basis_function[base1].gradient(point_2[l],gv_2),triangle_basis_function[base2].gradient(point_2[l],gv_2)));
			rhs(index4[base1])+=Jxy_2*f(point_2[l])*triangle_basis_function[base1].value(point_2[l],gv_2);
		}
		}
		/// TO DO: 计算每个积分点上的基函数梯度值，数值积分，拼装局部刚度矩阵，累加至整体刚度矩阵。
	}
	///处理所有边界条件
	for(unsigned int index=0;index < dim;index++)
	{
		if(nozeroperow[index]==3||nozeroperow[index]==4)
		{	
			int x_num=index%(n+1);
			int y_num=index/(n+1);
			double x=x_num*h;	
			double y=y_num*h;
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
	Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);	
	std::ofstream fs;
	///　输出到output.m用Matlab或Octave运行，得到计算结果
	fs.open("output.m");
	fs<<"x=0:1/"<<n<<":1;"<<std::endl;
	fs<<"y=0:1/"<<n<<":1;"<<std::endl;
	fs<<"[X,Y]=meshgrid(x,y);"<<std::endl;
	fs<<"U=[";
	int c=0;
	for(int j=0;j<n+1;j++)
	{	
		for(int i=0;i<n+1;i++)
			fs<<solution[i+c+n*j]<<" , ";
		fs<<";"<<std::endl;
		c++;
	}
	fs<<"]"<<std::endl;
	fs<<"surf(x,y,U);"<<std::endl;
    return 0;
};
