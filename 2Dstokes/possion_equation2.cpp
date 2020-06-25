/**
 * @file possion_equation.cpp
 * @author Lishijie (lsj1018845759@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2020-06-23
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <unordered_map>
#include <map>
#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/SparseMatrixTool.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/precondition.h>
#include <lac/solver_minres.h>
#include <lac/solver_gmres.h>
#include "../include/RectangleDomain.h"

#define PI (4.0*atan(1.0))
//定义边界条件；
double u(const double * p)
{
    return 1.0;
};

double f(const double * p)
{
    return 5 * PI * PI * u(p);
};

double zfun(const double * p)
{
    return 0.0;
};
int main(int argc, char* argv[])
{
	TemplateGeometry<2> rectangle_template_geometry2;
    rectangle_template_geometry2.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform2;
    rectangle_coord_transform2.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof2(rectangle_template_geometry2);
    rectangle_template_dof2.readData("rectangle.lagrange.2.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function2(rectangle_template_dof2);
    rectangle_basis_function2.readData("rectangle.lagrange.2.bas_fun");
    TemplateElement<double, 2, 2> template_element2;
    template_element2.reinit(rectangle_template_geometry2,
			    rectangle_template_dof2,
			    rectangle_coord_transform2,
			    rectangle_basis_function2);
	double volume2 = template_element2.volume();
    /// 取了 4 次代数精度。
	std::cout<<"Done!"<<std::endl;
    const QuadratureInfo<2>& quad_info2 = template_element2.findQuadratureInfo(4);   
    int n_quadrature_point2 = quad_info2.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point2 = quad_info2.quadraturePoint();
    int n_element_dof2 = template_element2.n_dof();
    int n_bas2 = rectangle_basis_function2.size();

    std::vector<AFEPack::Point<2> > gv2(9);
    std::vector<AFEPack::Point<2> > lv2(9);
	lv2[0][0] = -1.0;lv2[0][1] = -1.0;
	lv2[1][0] = 1.0;lv2[1][1] = -1.0;
	lv2[2][0] = 1.0;lv2[2][1] = 1.0;
	lv2[3][0] = -1.0;lv2[3][1] = 1.0;
	lv2[4][0] = 0.0;lv2[4][1] = -1.0;
	lv2[5][0] = 1.0;lv2[5][1] = 0.0;
	lv2[6][0] = 0.0;lv2[6][1] = 1.0;
	lv2[7][0] = -1.0;lv2[7][1] = 0.0;
	lv2[8][0] = 0.0;lv2[8][1] = 0.0;

    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform;
    rectangle_coord_transform.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof(rectangle_template_geometry);
    /// 一次元。
    rectangle_template_dof.readData("rectangle.1.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function(rectangle_template_dof);
    rectangle_basis_function.readData("rectangle.1.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    template_element.reinit(rectangle_template_geometry,
			    rectangle_template_dof,
			    rectangle_coord_transform,
			    rectangle_basis_function);
    /// 取了 4 次代数精度。
	double volume = template_element.volume();
    const QuadratureInfo<2>& quad_info = template_element.findQuadratureInfo(4);   
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point = quad_info.quadraturePoint();
    int n_element_dof = template_element.n_dof();
    int n_bas = rectangle_basis_function.size();
    std::vector<AFEPack::Point<2> > gv(4);
    std::vector<AFEPack::Point<2> > lv(4);
	lv[0][0] = -1.0;lv[0][1] = -1.0;
	lv[1][0] = 1.0;lv[1][1] = -1.0;
	lv[2][0] = 1.0;lv[2][1] = 1.0;
	lv[3][0] = -1.0;lv[3][1] = 1.0;

	
	/// 设置边界
    double x0 = 0.0;	
    double y0 = 0.0;
    double x1 = 1.0;
    double y1 = 1.0;
	int n = 20;
	int nx = n;
	int ny = n;
	//利用include里的函数RectangleDomain<2>对网格参数进行初始化；
	RectangleDomain<2> domain;
	/// 利用网格自动生成函数；
	domain.initial_2D_rectangle_domain(x0,x1,y0,y1,nx,ny);
	/// 生成网格，可以用domain.Mesh进行操作；由于处于测试版本，Mesh为public可改进；
	/// 这里网格生成的默认格式为四边形(二维情况)或六面体(三维情况)；可改进；
	domain.generate_mesh();
	/// 设置剖分断数和节点总数
	int dim = ( domain.get_Divide( 0 ) + 1 ) * ( domain.get_Divide( 1 ) + 1 );
    /// nozeroperow中每一个值表示对应行数非零元素个数
	RectangleDomain<2> domain2;
	domain2.initial_2D_rectangle_domain(x0,x1,y0,y1,nx,ny);
	domain2.set_divide_mode("Q2");
	domain2.generate_mesh();
	int dim2=( 2 * domain2.get_Divide( 0 ) + 1 ) * ( 2 * domain2.get_Divide( 1 ) + 1 );
	int dimension= 2*dim2+dim;
	std::vector<unsigned int> nozeroperow(dimension);
		for(int i = 0;i <2 * dim2;i++)
		nozeroperow[i]=15;
	nozeroperow[0]=nozeroperow[2 * n]=nozeroperow[2 * n * (2 * n + 1)]=nozeroperow[dim2 - 1] = 9;
	nozeroperow[dim2]=nozeroperow[dim2 + 2 * n]=nozeroperow[dim2 + 2 * n * (2 * n + 1)]=nozeroperow[2  *dim2 - 1] = 9;
	for(int i = 1;i < 2 * n;i = i + 2)
	{
		nozeroperow[i] = 9;
		nozeroperow[2 * n * (2 * n + 1) +  i] = 9;
		nozeroperow[dim2 + i]=9;
		nozeroperow[dim2 + 2 * n * (2 * n + 1) +  i] = 9;

	}
	for(int i = 1;i < 2 * n;i = i + 2)
	{
		nozeroperow[(2 * n + 1) * i] = 9;
		nozeroperow[(2 * n + 1) * i + 2 * n] = 9;
		nozeroperow[dim2 + (2 * n + 1) * i] = 9;
		nozeroperow[dim2 + (2 * n + 1) * i + 2 * n] = 9;
	}
	for(int i = 1 ;i < 2 * n;i = i + 2)
		for(int j = 1;j < 2 * n;j = j + 2)
		{
			nozeroperow[(2 * n + 1) * i + j] = 9;
			nozeroperow[dim2 + (2 * n + 1) * i + j] = 9;
		}
	for(int i = 2 ;i < 2 * n ;i = i + 2)
		for(int j = 2 ;j < 2 * n ;j = j + 2)
		{
			nozeroperow[i *(2 * n + 1) + j] = dim2;
			nozeroperow[dim2 + i *(2 * n + 1) + j] = dim2;
		}
	for(int i =2 * dim2;i < dimension;i++)
		nozeroperow[i]=2 * dim2;
	for(int i = 0;i < 2 * dim2 + dim;i++)
		nozeroperow[i]=nozeroperow[i]+dim;
		
	
	/// 请自行忽略这些奇怪的操作；
	SparsityPattern sp_stiff_matrix(dimension, nozeroperow);
	//使用SparsityPattern来定义SparseMatrix,在这里为了确定刚度矩阵非零元素的位置；

	for(element_iterator it_b=domain2.Mesh.begin();it_b!=domain2.Mesh.end();++it_b)
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
			for(index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
			{
				sp_stiff_matrix.add(it->second,it2->second);
				sp_stiff_matrix.add(it->second + dim2,it2->second + dim2);
			}
	element_iterator it_=domain.Mesh.begin();
	element_iterator it_2=domain2.Mesh.begin();
	element_iterator it_ee=domain.Mesh.end();
	element_iterator it_ee2=domain2.Mesh.end();
	while(it_!=it_ee||it_2!=it_ee2)
	{
		for(index_iterator it=it_2->begin();it!=it_2->end();++it)
			{
				for(index_iterator it2=it_->begin();it2!=it_->end();++it2)
				{	
					sp_stiff_matrix.add(it->second,it2->second + 2 * dim2);
					sp_stiff_matrix.add(it->second + dim2,it2->second + 2 * dim2);
				}
			}
			for(index_iterator it=it_->begin();it!=it_->end();++it)
			{
				for(index_iterator it2=it_2->begin();it2!=it_2->end();++it2)
				{
					sp_stiff_matrix.add(it->second + 2 * dim2,it2->second);
					sp_stiff_matrix.add(it->second + 2 * dim2,it2->second + dim2);
				}
			}
		++it_;
		++it_2;
	}	
	sp_stiff_matrix.compress();
	///下面组装刚度矩阵(assemble the stifff matrix)和右端项(rhs);
	SparseMatrix<double> stiff_mat(sp_stiff_matrix);
	//FullMatrix<double> stiff_mat(dimension,dimension);
	Vector<double> rhs( dimension );
	for(element_iterator it_b2=domain2.Mesh.begin();it_b2!=domain2.Mesh.end();++it_b2)
	{       
		for(index_iterator it=it_b2->begin();it!=it_b2->end();++it)
		{
			gv2[it->first.second][0]=it->first.first[0];
			gv2[it->first.second][1]=it->first.first[1];
		}
		for (int l = 0; l < n_quadrature_point2; l++)
		{
			std::vector<AFEPack::Point<2> > point2=rectangle_coord_transform2.local_to_global(q_point2, lv2, gv2);
			double Jxy=quad_info2.weight(l)*rectangle_coord_transform2.local_to_global_jacobian(q_point2[l], lv2, gv2)*volume2;
			for(index_iterator it=it_b2->begin();it!=it_b2->end();it++)
			{
				for(index_iterator it2=it_b2->begin();it2!=it_b2->end();it2++)
				{	
					stiff_mat.add(it->second,it2->second,Jxy*innerProduct(rectangle_basis_function2[it->first.second].gradient(point2[l],gv2),rectangle_basis_function2[it2->first.second].gradient(point2[l],gv2)));
					stiff_mat.add(it->second + dim2,it2->second + dim2,Jxy*innerProduct(rectangle_basis_function2[it->first.second].gradient(point2[l],gv2),rectangle_basis_function2[it2->first.second].gradient(point2[l],gv2)));
				}
				rhs(it->second)+=Jxy*f(point2[l])*rectangle_basis_function2[it->first.second].value(point2[l],gv2);
				rhs(it->second + dim2)+=Jxy*f(point2[l])*rectangle_basis_function2[it->first.second].value(point2[l],gv2);
			}
		}	
	}
	element_iterator it_b=domain.Mesh.begin();
	element_iterator it_b2=domain2.Mesh.begin();
	element_iterator it_e=domain.Mesh.end();
	element_iterator it_e2=domain2.Mesh.end();
	while(it_b!=it_e||it_b2!=it_e2)
	{       
		if(it_b!=it_e)
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
		{
			gv[it->first.second][0]=it->first.first[0];
			gv[it->first.second][1]=it->first.first[1];

		}
		if(it_b2!=it_e2)
		for(index_iterator it=it_b2->begin();it!=it_b2->end();++it)
		{
			gv2[it->first.second][0]=it->first.first[0];
			gv2[it->first.second][1]=it->first.first[1];

		}

		for (int l = 0; l < n_quadrature_point2; l++)
		{	
			std::vector<AFEPack::Point<2> > point=rectangle_coord_transform2.local_to_global(q_point, lv, gv);
			double Jxy=quad_info.weight(l)*rectangle_coord_transform2.local_to_global_jacobian(q_point[l], lv, gv)*volume;
			std::vector<AFEPack::Point<2> > point2=rectangle_coord_transform2.local_to_global(q_point2, lv2, gv2);
			double Jxy2=quad_info2.weight(l)*rectangle_coord_transform2.local_to_global_jacobian(q_point2[l], lv2, gv2)*volume2;

			for(index_iterator it=it_b2->begin();it!=it_b2->end();++it)
			{
				
				for(index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
				{	
					stiff_mat.add(it->second,it2->second + 2 * dim2,-Jxy*rectangle_basis_function2[it->first.second].value(point2[l],gv2)*rectangle_basis_function[it2->first.second].gradient(point[l],gv)[0]);
					stiff_mat.add(it->second + dim2,it2->second + 2 * dim2,-Jxy*rectangle_basis_function2[it->first.second].value(point2[l],gv2)*rectangle_basis_function[it2->first.second].gradient(point[l],gv)[1]);
				}
			}
			for(index_iterator it=it_b->begin();it!=it_b->end();++it)
			{
				for(index_iterator it2=it_b2->begin();it2!=it_b2->end();++it2)
				{
					stiff_mat.add(it->second + 2 * dim2,it2->second,-Jxy*rectangle_basis_function2[it2->first.second].value(point2[l],gv2)*rectangle_basis_function[it->first.second].gradient(point[l],gv)[0]);
					stiff_mat.add(it->second + 2 * dim2,it2->second + dim2,-Jxy*rectangle_basis_function2[it2->first.second].value(point2[l],gv2)*rectangle_basis_function[it->first.second].gradient(point[l],gv)[1]);
				}
			}
		}
		++it_b;
		++it_b2;	
	}
	
	for (int i = 0; i < dim2; i++)
    {	
    	if (nozeroperow[i] == 9 + dim||nozeroperow[i] == 15 + dim)
    	{
    	    SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(i);
    	    SparseMatrix<double>::iterator row_end = stiff_mat.end(i);
    	    double diag = row_iterator->value();
			double bnd_value = 0.0;
            rhs(i) = diag * bnd_value;
    	    for (++row_iterator; row_iterator != row_end; ++row_iterator)
            {
            	row_iterator->value() = 0.0;
    			int k = row_iterator->column();
                SparseMatrix<double>::iterator col_iterator = stiff_mat.begin(k);   
                SparseMatrix<double>::iterator col_end = stiff_mat.end(k);   
    	    	for (++col_iterator; col_iterator != col_end; ++col_iterator)
		    if (col_iterator->column() == i)
			break;
    		if (col_iterator == col_end)
    		{
		    std::cerr << "Error!" << std::endl;
		    exit(-1);
    		}
    		rhs(k) -= col_iterator->value() * bnd_value; 
    		col_iterator->value() = 0.0;	
            }

    	    row_iterator = stiff_mat.begin(i + dim2);
    	    row_end = stiff_mat.end(i + dim2);
    	    diag = row_iterator->value();
    	    //bnd_value = zfun(dof.interp_point); 
            rhs(i + dim2) = diag * bnd_value;
    	    for (++row_iterator; row_iterator != row_end; ++row_iterator)
            {
            	row_iterator->value() = 0.0;
    	    	int k = row_iterator->column();
                SparseMatrix<double>::iterator col_iterator = stiff_mat.begin(k);   
                SparseMatrix<double>::iterator col_end = stiff_mat.end(k);   
    	    	for (++col_iterator; col_iterator != col_end; ++col_iterator)
				{
		    		if (col_iterator->column() == i + dim2)
					break;
				}
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
	int n_v = dim2;
	int n_p = dim;
	int total_n_dof = dimension;
	std::vector<unsigned int> max_couple_vv(n_v);
    std::vector<unsigned int> max_couple_vp(n_v);
    
    for (int i = 0; i < n_v; i++)
    {
    	SparsityPattern::iterator row_start = sp_stiff_matrix.begin(i);
    	SparsityPattern::iterator row_end = sp_stiff_matrix.end(i);
    	for (; row_start != row_end; ++row_start)
    	{
	    int col = row_start->column();
    	    if (col < n_v)
		max_couple_vv[i]++;
	    else if (col >= 2 * n_v && col < total_n_dof)
		max_couple_vp[i]++;
     	}
    }
    SparsityPattern sp_vv(n_v, max_couple_vv);
    SparsityPattern sp_vp(n_v, n_p, max_couple_vp);
    for (int i = 0; i < n_v; i++)
    {
    	SparsityPattern::iterator row_start = sp_stiff_matrix.begin(i);
    	SparsityPattern::iterator row_end = sp_stiff_matrix.end(i);
	std::vector<unsigned int> vv_cols(max_couple_vv[i]);
	std::vector<unsigned int>::iterator vv_cols_iterator = vv_cols.begin();
	std::vector<unsigned int>::iterator vv_cols_start = vv_cols_iterator;
	std::vector<unsigned int> vp_cols(max_couple_vp[i]);
	std::vector<unsigned int>::iterator vp_cols_iterator = vp_cols.begin();
	std::vector<unsigned int>::iterator vp_cols_start = vp_cols_iterator;
    	for (; row_start != row_end; ++row_start)
    	{
	    int col = row_start->column();
    	    if (col < n_v)
	    {
		*vv_cols_iterator = col;
		++vv_cols_iterator;
	    }
	    else if (col >= 2 * n_v && col < total_n_dof)
	    {
		*vp_cols_iterator = col - 2 * n_v;
		++vp_cols_iterator;
	    }
	}
	sp_vv.add_entries(i, vv_cols_start, vv_cols_iterator);
	sp_vp.add_entries(i, vp_cols_start, vp_cols_iterator);
    }
    sp_vv.compress();
    sp_vp.compress();
    SparseMatrix<double> vxvx(sp_vv);
    SparseMatrix<double> vxp(sp_vp);
    SparseMatrix<double> vyvy(sp_vv);
    SparseMatrix<double> vyp(sp_vp);

    for (int i = 0; i < n_v; i++)
    {
    	SparseMatrix<double>::iterator row_start = stiff_mat.begin(i);
    	SparseMatrix<double>::iterator row_end = stiff_mat.end(i);
	std::vector<unsigned int> vv_cols(max_couple_vv[i]);
	std::vector<unsigned int>::iterator vv_cols_iterator = vv_cols.begin();
	std::vector<unsigned int>::iterator vv_cols_start = vv_cols_iterator;
	std::vector<double> vv_vals(max_couple_vv[i]);
	std::vector<double>::iterator vv_vals_iterator = vv_vals.begin();
	std::vector<double>::iterator vv_vals_start = vv_vals_iterator;
	std::vector<unsigned int> vp_cols(max_couple_vp[i]);
	std::vector<unsigned int>::iterator vp_cols_iterator = vp_cols.begin();
	std::vector<unsigned int>::iterator vp_cols_start = vp_cols_iterator;
	std::vector<double> vp_vals(max_couple_vp[i]);
	std::vector<double>::iterator vp_vals_iterator = vp_vals.begin();
	std::vector<double>::iterator vp_vals_start = vp_vals_iterator;
    	for (; row_start != row_end; ++row_start)
    	{
	    int col = row_start->column();
	    double val = row_start->value();
    	    if (col < n_v)
	    {
		*vv_cols_iterator = col;
		++vv_cols_iterator;
		*vv_vals_iterator = val;
		++vv_vals_iterator;
	    }
	    else if (col >= 2 * n_v && col < total_n_dof)
	    {
		*vp_cols_iterator = col - 2 * n_v;
		++vp_cols_iterator;
		*vp_vals_iterator = val;
		++vp_vals_iterator;

	    }
	}
	vxvx.add(i, vv_cols, vv_vals);
	vxp.add(i, vp_cols, vp_vals);
    }
    for (int i = n_v; i < 2 * n_v; i++)
    {
    	SparseMatrix<double>::iterator row_start = stiff_mat.begin(i);
    	SparseMatrix<double>::iterator row_end = stiff_mat.end(i);
	std::vector<unsigned int> vv_cols(max_couple_vv[i - n_v]);
	std::vector<unsigned int>::iterator vv_cols_iterator = vv_cols.begin();
	std::vector<unsigned int>::iterator vv_cols_start = vv_cols_iterator;
	std::vector<double> vv_vals(max_couple_vv[i - n_v]);
	std::vector<double>::iterator vv_vals_iterator = vv_vals.begin();
	std::vector<double>::iterator vv_vals_start = vv_vals_iterator;
	std::vector<unsigned int> vp_cols(max_couple_vp[i - n_v]);
	std::vector<unsigned int>::iterator vp_cols_iterator = vp_cols.begin();
	std::vector<unsigned int>::iterator vp_cols_start = vp_cols_iterator;
	std::vector<double> vp_vals(max_couple_vp[i - n_v]);
	std::vector<double>::iterator vp_vals_iterator = vp_vals.begin();
	std::vector<double>::iterator vp_vals_start = vp_vals_iterator;
    	for (; row_start != row_end; ++row_start)
    	{
	    int col = row_start->column();
	    double val = row_start->value();
    	    if (col >= n_v && col < 2 * n_v)
	    {
		*vv_cols_iterator = col - n_v;
		++vv_cols_iterator;
		*vv_vals_iterator = val;
		++vv_vals_iterator;
	    }
	    else if (col >= 2 * n_v && col < total_n_dof)
	    {
		*vp_cols_iterator = col - 2 * n_v;
		++vp_cols_iterator;
		*vp_vals_iterator = val;
		++vp_vals_iterator;
	    }
	}
	vyvy.add(i - n_v, vv_cols, vv_vals);
	vyp.add(i - n_v, vp_cols, vp_vals);
    }

	std::vector<unsigned int> max_couple_pp(n_p);
	for(element_iterator it_b=domain.Mesh.begin();it_b!=domain.Mesh.end();++it_b)
	{       
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
			max_couple_pp[it->second]+=n_element_dof;
	}
	SparsityPattern sp_pp(n_p, max_couple_pp);
	for (element_iterator it_b=domain.Mesh.begin();it_b!=domain.Mesh.end();++it_b) 
    {
	for (index_iterator it=it_b->begin();it!=it_b->end();++it)
	    for (index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
		sp_pp.add(it->second,it->second);
    }
	sp_pp.compress();
    SparseMatrix<double> mass_p(sp_pp);
	for(element_iterator it_b=domain.Mesh.begin();it_b!=domain.Mesh.end();++it_b) 
    {
		for(index_iterator it=it_b->begin();it!=it_b->end();++it)
		{
			gv[it->first.second][0]=it->first.first[0];
			gv[it->first.second][1]=it->first.first[1];
		}
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			std::vector<AFEPack::Point<2> > point=rectangle_coord_transform2.local_to_global(q_point, lv, gv);
			double Jxy=quad_info.weight(l)*rectangle_coord_transform2.local_to_global_jacobian(q_point[l], lv, gv)*volume;
			for(index_iterator it=it_b->begin();it!=it_b->end();++it)
				for(index_iterator it2=it_b->begin();it2!=it_b->end();++it2)
					mass_p.add(it->second, it2->second,Jxy*rectangle_basis_function[it->first.second].value(point[l],gv)*rectangle_basis_function[it2->first.second].value(point[l],gv));
		}
	}
	//for(int i = 0; i < dimension;i++)
		//std::cout<<rhs[i]<<std::endl;	

	//stiff_mat.print_formatted(std::cout,2,false);
	//AMGSolver solver(stiff_mat);
    /// 这里设置线性求解器的收敛判定为机器 epsilon 乘以矩阵的阶数，也
    /// 就是自由度总数。这个参数基本上是理论可以达到的极限。
    //solver.solve(solution, rhs, tol, 10000);
	/*std::ofstream fs;
	///　输出到output.m用Matlab或Octave运行，得到计算结果
	fs.open("output_u.m");
	fs<<"x="<<x0<<":1/"<<2 * n<<":"<<x1<<";"<<std::endl;
	fs<<"y="<<y0<<":1/"<<2 * n<<":"<<y1<<";"<<std::endl;
	//fs<<"[X,Y]=meshgrid(x,y);"<<std::endl;
	fs<<"X=[";
	for(int j=0;j<=2*n;j++)
	{	
		for(int i=0;i<=2*n;i++)
			fs<<solution[j * (2 * n + 1) + i]<<" , ";
		fs<<";"<<std::endl;
		;
	}
	fs<<"Y=[";
	for(int j=0;j<=2*n;j++)
	{	
		for(int i=0;i<=2*n;i++)
			fs<<solution[dim2 + j * (2 * n + 1) + i]<<" , ";
		fs<<";"<<std::endl;
		;
	}
	fs<<"]"<<std::endl;
	//fs<<"surf(x,y,U);"<<std::endl;*/
};
