使用步骤：
　　	矩形二次有限元重新定义了与本程序匹配的lagrange形基函数，
	定义在rectangle.lagrange.2文件夹中
	在运行程序之前要将rectangle.lagrange.2文件夹复制到/usr/local/AFEPack/template/rectangle中　
	然后在/usr/local/AFEPack/template/rectangle打开终端输入
	make
	即可．配置完成．
	
	
	随后在下载目录下打开regularMesh_Q2进行以下操作
	1:在终端打开.
	2:输入make clean
      	　　　make
      	　　./run
	3.在matlab中打开output.m运行得到计算结果
	说明：
  			1、利用std::pair<double,double>每一个节点的坐标；
  			   一维问题可以直接代替为double，三维问题目前初步考虑用std::vector<double>代替;
  			2、利用std::pair<std::pair<double,double>,unsigned int> 存储每一个节点坐标和局部编号；
  			3、利用std::map<std::pair<std::pair<double,double>,unsigned int>,unsigned int> 
  			   建立局部信息和整体编号之间的对应关系；
  			4、利用std::vector<std::map<std::pair<std::pair<double,double>,unsigned int>,unsigned int> >
  			   存储所有单元局部信息和整体编号之前的对应关系等信息；
  			5、定义宏函数local_coord_index_gen(n,x0,x1,y0,y1,i,j,local_index)，生成局部信息和整体编号之前的对应关系，
  			   其中n为网划分段数(可扩展为n1、n2、n3表示三维问题),x0 ~ y1为计算区域，i,j为单元里的位置参数；
  			   local_index为单元里局部编号；
  			6、定义宏函数global_to_coord(global_index,n,h)，计算整体编号和节点坐标之间的关系，其中global_index为整体编号；
  			7、定义宏函数可以减少结构化代码重复次数，同时可以提高程序的复用性和可读性；
  			8、定义宏函数之后Q1和Q2的实现方式可以趋于一致化，唯一不同的就是148-163行；
   			9、最后手动计算L2误差，by王老师。
