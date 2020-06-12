/**
 * @file   possion_equation.cpp
 * @author lsj <lsj@lsj>
 * @date   Tue Jun  12 19:05:19 2020
 * 
 * @brief  尝试将 AFEPack 采用线段线性有限元，混合边界条件，计算Possion方程
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


typedef std::unordered_map<unsigned int,int> index_map;

#define PI (4.0*atan(1.0))
double u(const double * p)
{
    //return sin(2 * PI * p[0]);
	//return exp( 2 * p[0] + p[1] );
	return 0;
};

double f(const double * p)
{
    return 2;
	//return -5 * u( p );
}; 
int main(int argc, char* argv[])
{
 
    TemplateGeometry<1> interval_template_geometry;
    interval_template_geometry.readData("interval.tmp_geo");
    CoordTransform<1, 1> interval_coord_transform;
    interval_coord_transform.readData("interval.crd_trs");
    TemplateDOF<1> interval_template_dof(interval_template_geometry);
    interval_template_dof.readData("interval.1.tmp_dof");
    BasisFunctionAdmin<double, 1, 1> interval_basis_function(interval_template_dof);
    interval_basis_function.readData("interval.1.bas_fun");
    TemplateElement<double, 1, 1> template_element;
    template_element.reinit(interval_template_geometry,
			    interval_template_dof,
			    interval_coord_transform,
			    interval_basis_function);

    double volume = template_element.volume();
    const QuadratureInfo<1>& quad_info = template_element.findQuadratureInfo(2);
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<1> > q_point = quad_info.quadraturePoint();
    int n_element_dof = template_element.n_dof();
    int n_bas = interval_basis_function.size();

	std::vector<AFEPack::Point<1> > arr(2);
    std::vector<AFEPack::Point<1> > gv(2);
    std::vector<AFEPack::Point<1> > lv(2);
	
	arr[0] = -1.0;
	arr[1] = 1.0;

    double xa = 0.5;	
    double xb = 1.5;
    int n = 10;
    double h = (xb - xa) / n;
	int dim = n + 1 ;
	Vector<double> rhs(dim);
	std::vector<unsigned int> nozeroperow(dim);
	for(int i =0;i<dim;i++)
		nozeroperow[i] = 3 ;
	nozeroperow[0]=nozeroperow[dim-1]=2;
	SparsityPattern sp_stiff_matrix(dim, nozeroperow);
	 for (int i = 0; i < n; i++)
    {
	    int idx0 = i ;
	    int idx1 = i + 1 ;
		index_map index({{0,idx0},{1,idx1}},n_element_dof);
		for(int base1 = 0;base1 < template_element.n_dof(); base1++)
			for(int base2 = 0;base2 < template_element.n_dof(); base2++)
				sp_stiff_matrix.add(index[base1],index[base2]);
	}
	sp_stiff_matrix.compress();
	SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    for (int i = 0; i < n; i++)
	{
	    double x0 = ((n - i) * xa + i * xb) / n;
	    int idx0 = i ;
	    double x1 = ((n - i - 1) * xa + ( i + 1 ) * xb) / n;
	    int idx1 = i + 1 ;
	    
	    int ele_idx =  i;

	    gv[0] = x0;
	    gv[1] = x1;
	    lv[0] = arr[0];
	    lv[1] = arr[1];
	    std::cout << ele_idx << ": " << std::endl;
	    std::cout << idx0 << ":(" << x0 << ") -> "
	    	      << idx1 << ":(" << x1 << ") "<< std::endl;
	    /// 现在尝试输出具体每个单元的积分点。
		index_map index({{0,idx0},{1,idx1}},n_element_dof);
	    for (int l = 0; l < n_quadrature_point; l++)
		{
		auto point=interval_coord_transform.local_to_global(q_point, lv, gv);
		double Jxy=quad_info.weight(l)*interval_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)*volume;
		for(int base1 = 0;base1 < template_element.n_dof(); base1++)
		{
			for(int base2 = 0;base2 < template_element.n_dof(); base2++)
				stiff_mat.add(index[base1],index[base2],Jxy*innerProduct(interval_basis_function[base1].gradient(point[l],gv),interval_basis_function[base2].gradient(point[l],gv)));
			rhs(index[base1])+=Jxy*f(point[l])*interval_basis_function[base1].value(point[l],gv);
		}
		}
	}
	for(unsigned int index=0;index < dim;index++)
	{
		if(index==0)
		{	
			double x = ((n - index) * xa + index * xb) / n;
			SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(index);
    	    SparseMatrix<double>::iterator row_end = stiff_mat.end(index);
    	    double diag = row_iterator->value();
			AFEPack::Point<1> bnd_point;
			bnd_point[0]=x;
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
	AMGSolver solver(stiff_mat);
	Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);	
	std::ofstream fs;
	///　输出到output.m用Matlab或Octave运行，得到计算结果
	fs.open("output.m");
	fs<<"x="<<xa<<"0:1/"<<n<<":"<<xb<<";"<<std::endl;
	fs<<"U=[";
	for(int j=0;j<n+1;j++)
		fs<<solution[j]<<" , ";
	fs<<"]"<<std::endl;
	fs<<"plot(x,U);"<<std::endl;
    return 0;
	
};
