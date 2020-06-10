/**
 * @file   step2.cpp
 * @author Lsj <Lsj@Lsj>
 * @date   Tue Jun  5 9:12:34 2020
 * 
 * @brief  利用AFEPack里的函数以及功能，划分结构化网格，使用lagrange九节点二次基函数
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
    //return sin(4 * PI * p[0]) * sin(2 * PI * p[1]);
	return sin( 2 * PI * p[0]) * sin(PI * p[1]);
};

double f(const double * p)
{
    //return 20 * PI * PI * u(p);
	return 5 * PI * PI* u( p );
}; 

typedef std::unordered_map<unsigned int,int> index_map;

int main(int argc, char* argv[])
{
    /// 这里基本上和 possion_equation 中配置一致。对比
    /// possion_equation_manual 看更清楚。
    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangleregular.tmp_geo");
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

    /// 产生一个具体单元顶点的缓存。
    double ** arr = (double **) new double* [9];
    for (int i = 0; i < 9; i++)
	arr[i] = (double *) new double [2];
    std::vector<AFEPack::Point<2> > gv(9);
    std::vector<AFEPack::Point<2> > lv(9);

    for (int i = 0; i < n_element_dof; i++)
    {
	AFEPack::Point<2> pnt = q_point[i];
	arr[0][0] = -1.0;arr[0][1] = -1.0;
	arr[1][0] = 1.0;arr[1][1] = -1.0;
	arr[2][0] = 1.0;arr[2][1] = 1.0;
	arr[3][0] = -1.0;arr[3][1] = 1.0; 
	arr[4][0] = 0.0;arr[4][1] = -1.0;
	arr[5][0] = 1.0;arr[5][1] = 0.0;
	arr[6][0] = 0.0;arr[6][1] = 1.0;
	arr[7][0] = -1.0;arr[7][1] = 0.0;
	arr[8][0] = 0.0;arr[8][1] =0.0;

	lv[0][0] = arr[0][0];lv[0][1] = arr[0][1];
	lv[1][0] = arr[1][0];lv[1][1] = arr[1][1];
	lv[2][0] = arr[2][0];lv[2][1] = arr[2][1];
	lv[3][0] = arr[3][0];lv[3][1] = arr[3][1];
	lv[4][0] = arr[4][0];lv[4][1] = arr[4][1];
	lv[5][0] = arr[5][0];lv[5][1] = arr[5][1];
	lv[6][0] = arr[6][0];lv[6][1] = arr[6][1];
	lv[7][0] = arr[7][0];lv[7][1] = arr[7][1];
	lv[8][0] = arr[8][0];lv[8][1] = arr[8][1];
    }
    double xa = 0.0;	
    double ya = 0.0;
    double xb = 1.0;
    double yb = 1.0;
    int n = 10;
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
			index_map index({{0,idx0},{1,idx1},{2,idx2},{3,idx3},{4,idx4},{5,idx5},{6,idx6},{7,idx7},{8,idx8}},template_element.n_dof());
			for(int i = 0;i < template_element.n_dof(); i++)
				for(int j = 0;j < template_element.n_dof(); j++)
					sp_stiff_matrix.add(index[i],index[j]);
		}

	sp_stiff_matrix.compress();
	
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);

    double h = (xb - xa) / n;
    for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
	    	double x0 = ((2 * n - 2 * i) * xa + 2 * i * xb) / (2 * n);
	    	double y0 = ((2 * n - 2 * j) * ya + 2 * j * yb) / (2 * n);
	    	int idx0 = (2 * j) * (2 * n + 1) + (2 * i) ; 

			double x1 = ((2 * n  - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y1 = ((2 * n  - 2 * j) * ya + 2 * j * yb) / (2 * n);
	    	int idx1 = (2 * j) * (2 * n + 1) + (2 * i + 2) ;

	    	double x2 = ((2 * n - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y2 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
	    	int idx2 = (2 * j + 2) * (2 * n + 1) + (2 * i + 2) ;

	    	double x3 = ((2 * n - 2 * i) * xa + 2 * i * xb) / (2 * n);
	    	double y3 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
	   		int idx3 = (2 * j + 2) * (2 * n + 1) + (2 * i) ; 

            double x4 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y4 = ((2 * n - 2 * j) * ya + (2 * j ) * yb) / (2 * n);
	    	int idx4 = (2 * j) * (2 * n + 1) + (2 * i + 1) ;

			double x5 = ((2 * n - 2 * i - 2) * xa + (2 * i + 2) * xb) / (2 * n);
	    	double y5 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx5 = (2 * j + 1) * (2 * n + 1) + (2 * i + 2) ;

			double x6 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y6 = ((2 * n - 2 * j - 2) * ya + (2 * j + 2) * yb) / (2 * n);
			int idx6 = (2 * j + 2) * (2 * n + 1) + (2 * i + 1) ;

			double x7 = ((2 * n - 2 * i) * xa + (2 * i) * xb) / (2 * n);
	    	double y7 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx7 = (2 * j + 1) * (2 * n + 1) + (2 * i) ;

			double x8 = ((2 * n - 2 * i - 1) * xa + (2 * i + 1) * xb) / (2 * n);
	    	double y8 = ((2 * n - 2 * j - 1) * ya + (2 * j + 1) * yb) / (2 * n);
			int idx8 = (2 * j + 1) * (2 * n + 1) + (2 * i + 1) ;

	    gv[0][0] = x0;gv[0][1] = y0;  
	    gv[1][0] = x1;gv[1][1] = y1;
	    gv[2][0] = x2;gv[2][1] = y2;
	    gv[3][0] = x3;gv[3][1] = y3;
	    gv[4][0] = x4;gv[4][1] = y4;
	    gv[5][0] = x5;gv[5][1] = y5;
	    gv[6][0] = x6;gv[6][1] = y6;
	    gv[7][0] = x7;gv[7][1] = y7;
	    gv[8][0] = x8;gv[8][1] = y8;

	    // 现在尝试输出具体每个单元的积分点。
	   for (int l = 0; l < n_quadrature_point; l++)
		{
		auto point=rectangle_coord_transform.local_to_global(q_point, lv, gv);
		double Jxy=quad_info.weight(l)*rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)*volume;
		index_map index({{0,idx0},{1,idx1},{2,idx2},{3,idx3},{4,idx4},{5,idx5},{6,idx6},{7,idx7},{8,idx8}},template_element.n_dof());
		for(int base1 = 0;base1 < template_element.n_dof(); base1++)
		{
			for(int base2 = 0;base2 < template_element.n_dof(); base2++)
				stiff_mat.add(index[base1],index[base2],Jxy*innerProduct(rectangle_basis_function[base1].gradient(point[l],gv),rectangle_basis_function[base2].gradient(point[l],gv)));
			rhs(index[base1])+=Jxy*f(point[l])*rectangle_basis_function[base1].value(point[l],gv);
		}
	   
		}
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
			std::cout<<"("<<x<<","<<y<<")"<<std::endl;
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
};
