/**
 * @file   possion_equation.cpp
 * @author Li ShiJie <lsj@lsj>
 * @date   Tue Jun  13 18:11:56 2020
 * 
 * @brief  利用AFEPack里的函数以及功能，划分结构化网格，使用双线性四节点基函数
 * 		　　，实现Possion方程编值问题计算．
 * 			说明：
 * 			1、利用std::pair<double,double>每一个节点的坐标；
 * 			   一维问题可以直接代替为double，三维问题目前初步考虑用std::vector<double>代替;
 * 			2、利用std::pair<std::pair<double,double>,unsigned int> 存储每一个节点坐标和局部编号；
 * 			3、利用std::map<std::pair<std::pair<double,double>,unsigned int>,unsigned int> 
 * 			   建立局部信息和整体编号之间的对应关系；
 * 			4、利用std::vector<std::map<std::pair<std::pair<double,double>,unsigned int>,unsigned int> >
 * 			   存储所有单元局部信息和整体编号之前的对应关系等信息；
 * 			5、定义宏函数local_coord_index_gen(n,x0,x1,y0,y1,i,j,local_index)，生成局部信息和整体编号之前的对应关系，
 * 			   其中n为网划分段数(可扩展为n1、n2、n3表示三维问题),x0 ~ y1为计算区域，i,j为单元里的位置参数；
 * 			   local_index为单元里局部编号；
 * 			6、定义宏函数global_to_coord(global_index,n,h)，计算整体编号和节点坐标之间的关系，其中global_index为整体编号；
 * 			7、定义宏函数可以减少结构化代码重复次数，同时可以提高程序的复用性和可读性；
 *			8、最后手动计算L2误差，by王老师。
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
#include <AFEPack/EasyMesh.h>
#include <AFEPack/SparseMatrixTool.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>

#define PI (4.0*atan(1.0))
//定义边界条件；
double u(const double * p)
{
    return sin(PI * p[0]) * sin(2 * PI * p[1]);
};
//定义非齐次项函数；
double f(const double * p)
{
    return 5 * PI * PI * u(p);
}; 

typedef std::pair<double,double> coord_pair;
typedef std::pair<coord_pair,unsigned int > coord_loacl_pair;
typedef std::map<coord_loacl_pair,unsigned int> loacl_to_global_map;
typedef std::vector<loacl_to_global_map>  loacl_to_global_vector;

#define local_coord_index_gen(n,x0,x1,y0,y1,i,j,local_index){\
		double x = ((n - (i)) * x0 + i * x1) / n;\
	    double y = ((n - (j)) * y0 + j * y1) / n;\
		coord_pair coord(x,y);\
		coord_loacl_pair coord_loacl(coord,local_index);\
		loacl_global_map[coord_loacl]=(j) * (n + 1) + (i);\
		}

#define global_to_coord(global_index,n,h){\
		double x_coord=( index % ( n + 1 ) ) * h;\
		double y_coord=( index / ( n + 1 ) ) * h;\
		}

int main(int argc, char* argv[])
{
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

    double volume = template_element.volume();
    /// 取了 4 次代数精度。
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
	/// 设置剖分断数和节点总数
    int n = 100;
	double h = (x1 - x0) / n;
	int dim = ( n + 1 ) * ( n + 1 );
	Vector<double> rhs( dim );
    /// nozeroperow中每一个值表示对应行数非零元素个数
	std::vector<unsigned int> nozeroperow(dim);
	for(int i = 0;i <= dim ; i++)
		nozeroperow[ i ] = 9;
	nozeroperow[ 0 ] = 4;
	nozeroperow[ dim - 1 ] = 4;
	nozeroperow[ n ] = 4;
	nozeroperow[ dim - 1 - n] = 4;
	for(int i = 1;i < n;i++)
	{
		nozeroperow[ i ]=6;
		nozeroperow[ dim - 1 - i]=6;
		nozeroperow[ i * ( n + 1 ) ]=6;
		nozeroperow[( i + 1 ) * n + i ] = 6;
	}
	/// 填充非零元素对应的行索引和列索引
	loacl_to_global_vector loacl_global_vector(n * n);
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
		/// 现在尝试输出具体每个单元的积分点。
		///　合成整体刚度矩阵
		//* 6----7----8
		///	|	 |    |
		///	3----4----5     
    	//  |    |    |
		//	0----1----2
		//
		///element 1:0->1->4->3
		///element 2:1->2->5->4
		///element 3:3->4->7->6
		///element 4:4->5->8->7
		int ele_idx = j * n + i;
		loacl_to_global_map loacl_global_map;
		//此处利用之前所定义的宏local_coord_index_gen，生成每一个单元element对应的
		//局部坐标、编号和整体编号之前的对应关系，用std::map进行定义；
		local_coord_index_gen(n,x0,x1,y0,y1,i,j,0);
		local_coord_index_gen(n,x0,x1,y0,y1,i + 1,j,1);
		local_coord_index_gen(n,x0,x1,y0,y1,i + 1,j + 1,2);
		local_coord_index_gen(n,x0,x1,y0,y1,i,j + 1,3);
		//每一个std::map可以存储到std::vector中去，未尝试push_back（）函数；
		loacl_global_vector[ele_idx]=loacl_global_map;
	}
	SparsityPattern sp_stiff_matrix(dim, nozeroperow);
	//使用SparsityPattern来定义SparseMatrix,在这里为了确定刚度矩阵非零元素的位置；
	for(auto it_b=loacl_global_vector.begin();it_b!=loacl_global_vector.end();++it_b)
		for(auto it=(*it_b).begin();it!=(*it_b).end();++it)
			for(auto it2=(*it_b).begin();it2!=(*it_b).end();++it2)
				sp_stiff_matrix.add(it->second,it2->second);
	sp_stiff_matrix.compress();
	///下面组装刚度矩阵(assemble the stifff matrix)和右端项(rhs);
	SparseMatrix<double> stiff_mat(sp_stiff_matrix);
	for(auto it_b=loacl_global_vector.begin();it_b!=loacl_global_vector.end();++it_b)
	{       
		for(auto it=it_b->begin();it!=it_b->end();++it)
		{
			gv[it->first.second][0]=it->first.first.first;
			gv[it->first.second][1]=it->first.first.second;
		}
		for (int l = 0; l < n_quadrature_point; l++)
		{
			auto point=rectangle_coord_transform.local_to_global(q_point, lv, gv);
			double Jxy=quad_info.weight(l)*rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv)*volume;
			for(auto it=(*it_b).begin();it!=(*it_b).end();it++)
			{
				for(auto it2=(*it_b).begin();it2!=(*it_b).end();it2++)	
					stiff_mat.add(it->second,it2->second,Jxy*innerProduct(rectangle_basis_function[it->first.second].gradient(point[l],gv),rectangle_basis_function[it2->first.second].gradient(point[l],gv)));
				rhs(it->second)+=Jxy*f(point[l])*rectangle_basis_function[it->first.second].value(point[l],gv);
			}
		}	
	}
	///处理所有边界条件
	for(unsigned int index=0;index < dim;index++)
	{
		if(nozeroperow[index]==4||nozeroperow[index]==6)
		{	
			double x_coord,y_coord=0.0;
			global_to_coord(index,n,h);
			SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(index);
    	    SparseMatrix<double>::iterator row_end = stiff_mat.end(index);
    	    double diag = row_iterator->value();
			AFEPack::Point<2> bnd_point;
			bnd_point[0]=x_coord;
			bnd_point[1]=y_coord;
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
	fs<<"x="<<x0<<":1/"<<n<<":"<<x1<<";"<<std::endl;
	fs<<"y="<<y0<<":1/"<<n<<":"<<y1<<";"<<std::endl;
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

	double error = 0.0;
    for(unsigned int index = 0; index < dim; index++)
    {
	int x_num = index % (n + 1);
	int y_num = index / (n + 1);
	double x = x_num * h;	
	double y = y_num * h;
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
