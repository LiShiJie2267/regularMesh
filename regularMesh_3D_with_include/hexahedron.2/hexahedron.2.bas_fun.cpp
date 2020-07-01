/*************************************************************************
 * 
 */

#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

#define GAUSS_ELIMINATION 								\
	for (i = 0;i < 7;i ++) {							\
		k = i;									\
		for (j = i + 1;j < 8;j ++)						\
			if (fabs(m[j][i]) > fabs(m[k][i])) k = j;			\
		if (k != i) {								\
			for (j = i;j < 8;j ++) {					\
				tmp = m[i][j];						\
				m[i][j] = m[k][j];					\
				m[k][j] = tmp;						\
			}								\
			tmp = a[i][0];							\
			a[i][0] = a[k][0];						\
			a[k][0] = tmp;							\
			tmp = a[i][1];							\
			a[i][1] = a[k][1];						\
			a[k][1] = tmp;              \
			tmp = a[i][2];							\
			a[i][2] = a[k][2];						\
			a[k][2] = tmp;							\
		}									\
		for (j = i + 1;j < 8;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i + 1;k < 8;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];           \
			a[j][2] -= tmp*a[i][2];						\
		}									\
	}										\
	a[7][0] /= m[7][7];								\
	a[7][1] /= m[7][7];               \
	a[7][2] /= m[7][7];								\
	for (i = 6;i >= 0;i--) {							\
		for (j = i + 1;j < 8;j++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];          \
      a[i][2] -= m[i][j]*a[j][2];					\
		}									\
		a[i][0] /= m[i][i];							\  
		a[i][1] /= m[i][i];             \
    a[i][2] /= m[i][i];						\
	}

#define COMMON_PART                                                                     \
    int i, j, k;                                                                        \
    double m[8][8], tmp, xi, eta,theta;                                                       \
    double a[8][3] = {{-1.0, -1.0, -1.0},{1.0, -1.0,-1.0},{1.0, 1.0,-1.0},{-1.0, 1.0,-1.0},{-1.0, -1.0, 1.0},{1.0, -1.0,1.0},{1.0,1.0,1.0},{-1.0,1.0,1.0}}; \
    for (i = 0;i < 8;i ++) {                                                            \
      m[i][0] = 1.0;                                                                    \
      m[i][1] = v[i][0];                                                                \
      m[i][2] = v[i][1];                                                                \
      m[i][3] = v[i][2];                                                                \
      m[i][4] = v[i][0]*v[i][1];                                                        \
      m[i][5] = v[i][0]*v[i][2];                                                        \
      m[i][6] = v[i][1]*v[i][2];                                                        \
      m[i][7] = v[i][0]*v[i][1]*v[i][2];                                                \
    }                                                                                   \
    GAUSS_ELIMINATION;                                                                  \
    xi  = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[2]+ a[4][0]*p[0]*p[1] + a[5][0]*p[0]*p[2] + a[6][0]*p[1]*p[2] + a[7][0]*p[0]*p[1]*p[2];   \        
    eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[2]+ a[4][1]*p[0]*p[1] + a[5][1]*p[0]*p[2] + a[6][1]*p[1]*p[2] + a[7][1]*p[0]*p[1]*p[2];   \  
    theta = a[0][2] + a[1][2]*p[0] + a[2][2]*p[1] + a[3][2]*p[2]+ a[4][2]*p[0]*p[1] + a[5][2]*p[0]*p[2] + a[6][2]*p[1]*p[2] + a[7][2]*p[0]*p[1]*p[2];  \           

#define D_XI_DX  ( a[1][0] + a[4][0]*p[1] + a[5][0]*p[2] + a[7][0]*p[1]*p[2])
#define D_XI_DY  ( a[2][0] + a[4][0]*p[0] + a[6][0]*p[2] + a[7][0]*p[0]*p[2])
#define D_XI_DZ  ( a[3][0] + a[5][0]*p[0] + a[6][0]*p[1] + a[7][0]*p[0]*p[1])
#define D_ETA_DX  ( a[1][1] + a[4][1]*p[1] + a[5][1]*p[2] + a[7][1]*p[1]*p[2])
#define D_ETA_DY  ( a[2][1] + a[4][1]*p[0] + a[6][1]*p[2] + a[7][1]*p[0]*p[2])
#define D_ETA_DZ  ( a[3][1] + a[5][1]*p[0] + a[6][1]*p[1] + a[7][1]*p[0]*p[1])
#define D_THETA_DX  ( a[1][2] + a[4][2]*p[1] + a[5][2]*p[2] + a[7][2]*p[1]*p[2])
#define D_THETA_DY  ( a[2][2] + a[4][2]*p[0] + a[6][2]*p[2] + a[7][2]*p[0]*p[2])
#define D_THETA_DZ  ( a[3][2] + a[5][2]*p[0] + a[6][2]*p[1] + a[7][2]*p[0]*p[1])

  void lambda_1(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 - eta)*(1.0 - theta);
	  val[0]/=8;
  }

  void gradient_lambda_1(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (eta - 1.0)*(1.0 - theta) * D_XI_DX + (xi - 1.0)*(1.0 - theta) * D_ETA_DX + (xi - 1.0)*(1.0 - eta) * D_THETA_DX;
    val[1] = (eta - 1.0)*(1.0 - theta) * D_XI_DY + (xi - 1.0)*(1.0 - theta) * D_ETA_DY + (xi - 1.0)*(1.0 - eta) * D_THETA_DY;
    val[2] = (eta - 1.0)*(1.0 - theta) * D_XI_DZ + (xi - 1.0)*(1.0 - theta) * D_ETA_DZ + (xi - 1.0)*(1.0 - eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_2(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + xi)*(1.0 - eta)*(1.0 - theta);
	  val[0]/=8;
  }

  void gradient_lambda_2(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - eta)*(1.0 - theta) * D_XI_DX - (xi + 1.0)*(1.0 - theta) * D_ETA_DX - (xi + 1.0)*(1.0 - eta) * D_THETA_DX;
    val[1] = (1.0 - eta)*(1.0 - theta) * D_XI_DY - (xi + 1.0)*(1.0 - theta) * D_ETA_DY - (xi + 1.0)*(1.0 - eta) * D_THETA_DY;
    val[2] = (1.0 - eta)*(1.0 - theta) * D_XI_DZ - (xi + 1.0)*(1.0 - theta) * D_ETA_DZ - (xi + 1.0)*(1.0 - eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_3(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + xi)*(1.0 + eta)*(1.0 - theta);
    val[0]/=8;
  }

  void gradient_lambda_3(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + eta)*(1.0 - theta) * D_XI_DX + (xi + 1.0)*(1.0 - theta) * D_ETA_DX - (xi + 1.0)*(1.0 + eta) * D_THETA_DX;
    val[1] = (1.0 + eta)*(1.0 - theta) * D_XI_DY + (xi + 1.0)*(1.0 - theta) * D_ETA_DY - (xi + 1.0)*(1.0 + eta) * D_THETA_DY;
    val[2] = (1.0 + eta)*(1.0 - theta) * D_XI_DZ + (xi + 1.0)*(1.0 - theta) * D_ETA_DZ - (xi + 1.0)*(1.0 + eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_4(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 + eta)*(1.0 - theta);
    val[0]/=8;
  }

  void gradient_lambda_4(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = -(1.0 + eta)*(1.0 - theta) * D_XI_DX + (1.0 - xi)*(1.0 - theta) * D_ETA_DX - (1.0 - xi)*(1.0 + eta) * D_THETA_DX;
    val[1] = -(1.0 + eta)*(1.0 - theta) * D_XI_DY + (1.0 - xi)*(1.0 - theta) * D_ETA_DY - (1.0 - xi)*(1.0 + eta) * D_THETA_DY;
    val[2] = -(1.0 + eta)*(1.0 - theta) * D_XI_DZ + (1.0 - xi)*(1.0 - theta) * D_ETA_DZ - (1.0 - xi)*(1.0 + eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_5(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 - eta)*(1.0 + theta);
    val[0]/=8;
  }

  void gradient_lambda_5(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = -(1.0 - eta)*(1.0 + theta) * D_XI_DX - (1.0 - xi)*(1.0 + theta) * D_ETA_DX + (1.0 - xi)*(1.0 - eta) * D_THETA_DX;
    val[1] = -(1.0 - eta)*(1.0 + theta) * D_XI_DY - (1.0 - xi)*(1.0 + theta) * D_ETA_DY + (1.0 - xi)*(1.0 - eta) * D_THETA_DY;
    val[2] = -(1.0 - eta)*(1.0 + theta) * D_XI_DZ - (1.0 - xi)*(1.0 + theta) * D_ETA_DZ + (1.0 - xi)*(1.0 - eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_6(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + xi)*(1.0 - eta)*(1.0 + theta);
    val[0]/=8;
  }

  void gradient_lambda_6(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - eta)*(1.0 + theta) * D_XI_DX - (1.0 + xi)*(1.0 + theta) * D_ETA_DX + (1.0 + xi)*(1.0 - eta) * D_THETA_DX;
    val[1] = (1.0 - eta)*(1.0 + theta) * D_XI_DY - (1.0 + xi)*(1.0 + theta) * D_ETA_DY + (1.0 + xi)*(1.0 - eta) * D_THETA_DY;
    val[2] = (1.0 - eta)*(1.0 + theta) * D_XI_DZ - (1.0 + xi)*(1.0 + theta) * D_ETA_DZ + (1.0 + xi)*(1.0 - eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_7(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + xi)*(1.0 + eta)*(1.0 + theta);
    val[0]/=8;
  }

  void gradient_lambda_7(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 + eta)*(1.0 + theta) * D_XI_DX + (1.0 + xi)*(1.0 + theta) * D_ETA_DX + (1.0 + xi)*(1.0 + eta) * D_THETA_DX;
    val[1] = (1.0 + eta)*(1.0 + theta) * D_XI_DY + (1.0 + xi)*(1.0 + theta) * D_ETA_DY + (1.0 + xi)*(1.0 + eta) * D_THETA_DY;
    val[2] = (1.0 + eta)*(1.0 + theta) * D_XI_DZ + (1.0 + xi)*(1.0 + theta) * D_ETA_DZ + (1.0 + xi)*(1.0 + eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }

  void lambda_8(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 + eta)*(1.0 + theta);
    val[0]/=8;
  }

  void gradient_lambda_8(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = -(1.0 + eta)*(1.0 + theta) * D_XI_DX + (1.0 - xi)*(1.0 + theta) * D_ETA_DX + (1.0 - xi)*(1.0 + eta) * D_THETA_DX;
    val[1] = -(1.0 + eta)*(1.0 + theta) * D_XI_DY + (1.0 - xi)*(1.0 + theta) * D_ETA_DY + (1.0 - xi)*(1.0 + eta) * D_THETA_DY;
    val[2] = -(1.0 + eta)*(1.0 + theta) * D_XI_DZ + (1.0 - xi)*(1.0 + theta) * D_ETA_DZ + (1.0 - xi)*(1.0 + eta) * D_THETA_DZ;
    val[0]/=8;
    val[1]/=8;
    val[2]/=8;
  }
#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
