#include "Savitzky_Golay_2d.h"
#include <iostream>

using namespace std;

int Filter2d::alloc()
{
    filter =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    data = (double *) fftw_malloc(sizeof(double)*(size));
    datat = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    datat_filtered =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    plan_data2filter = fftw_plan_dft_r2c_2d(y_size, x_size, data, filter, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_data2datat =  fftw_plan_dft_r2c_2d(y_size, x_size, data, datat, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_filtered2data = fftw_plan_dft_c2r_2d(y_size, x_size, datat_filtered, data,FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    
}

Filter2d::Filter2d(int nx, int ny, int ld, int m, int xsize, int ysize):x_size(xsize),y_size(ysize)
{
//    double *temp = reinterpret_cast<double *>( fftw_malloc(sizeof(double) * size));
    int m1 = m+1;
    int x_len = nx*2 + 1;
    int y_len = ny*2 + 1;
	int col_len = x_len * y_len;
    size = x_size * y_size;
    size2 = (x_size/2 +1) * y_size;
    alloc();
    box_range.init(x_len, y_len);
    filter_range.init(xsize, ysize);


    int m11 = m1*(m1+1)/2;
    double *jacobian = new double[ m11 * col_len ];
    double *normal_matrix = new double[ m11 * m11];
	for(int k = -ny; k <= ny; k++){
        for(int s = -nx; s <= nx; s++){
            int ii =box_range(s, k);
            int col = 0;
            for (int i = 0; i < m1; i++){
                for(int j = 0; j < m1-i; j++, col++){
 //                   int jj =((i * m1) + j);
                    double temp = jacobian[ ii + col* col_len]
                    = pow (k, i) * pow(s, j);
                    //cout << ii <<" "<< col <<" "<<temp<<endl;
                }
            }
        }
        //cout << endl;
    }
    
	int mat_size = col_len*col_len;	
//	double *unit_mat = new double[ mat_size ];
//	for(int i = 0; i < mat_size; i+= col_len)  unit_mat[i] = 1;

	
	double double_one = 1.;
	double beta = 0.;
    int row_len = m11;
//	cout <<"dgemm"<<endl;
//	cout << col_len <<" "<<col_len<<" "<<row_len<<" "<<double_one<<" "<<beta<<endl;
	cblas_dgemm(CblasColMajor,  CblasTrans,  CblasNoTrans, row_len, row_len, col_len, double_one, jacobian, col_len, jacobian, col_len,
				beta, normal_matrix, row_len); // Normal = J*J
	//cout <<"dgemm finished"<<endl;
	
	int *ipiv = new int[row_len]; 
	int info;
    for(int i = 0; i < row_len; i++){
        for(int j = 0; j< row_len; j++){
            //cout << normal_matrix[j*row_len + i]<<" ";
        }
//        cout << endl;
    }
	
	
	dgetrf_(&row_len, &row_len, normal_matrix, &row_len, ipiv, &info); //rf normal
#ifdef VERBOSE_
    for(int i = 0; i < row_len; i++){
        for(int j = 0; j< row_len; j++){
            cout << normal_matrix[j*row_len + i]<<" ";
        }
        cout<< ipiv[i] << endl;
    }
    cout<<endl;
#endif
	double *work;
	double worksize;
	int lwork = -1;	
	dgetri_(&row_len, normal_matrix,&row_len, ipiv, &worksize, &lwork, &info);
	lwork = (int) (worksize);
	work = new double [lwork];
	
	dgetri_(&row_len, normal_matrix,&row_len, ipiv, work, &lwork, &info);
#ifdef VERBOSE_	
    for(int i = 0; i < row_len; i++){
        for(int j = 0; j< row_len; j++){
            cout << normal_matrix[j*row_len + i]<<" ";
        }
        cout<< ipiv[i] << endl;
    }
	cout << lwork << "dgetri"<<endl;
    for(int i = 0; i < col_len; i++){
        for(int j = 0; j < row_len; j++ ){
            cout << jacobian[j*col_len + i]<<" ";
        }
        cout << endl;
    }
#endif
	double *cof=new double [row_len * col_len];
	cblas_dgemm(CblasColMajor,  CblasNoTrans, CblasTrans, row_len, col_len, row_len, double_one, normal_matrix, row_len, jacobian, col_len,
				beta, cof, row_len);
    
#ifdef VERBOSE_
    cout <<"cof"<<endl;
    int posi = 0;
    for(int j = 0; j < y_len; j++){
        for(int i = 0; i < x_len; i++, posi++){
            cout << cof[posi*row_len]<<" ";
        }cout << endl;
    }
	cout <<"dgemm"<<endl;
	cout << "reset"<<endl;
#endif
    for(int i = 0; i < size; i++) data[i] = 0;
//	cout << "copie"<<endl;
    for(int j = -ny; j <= ny; j++){
        for(int i = -nx; i <= nx; i++) {
            double temp = data[filter_range(i, j)] //box_range?
            = cof[row_len*(box_range(i, j))+ld];
    //        cout << temp<<" "<<box_range(i, j)<<" ";
        }
        //cout << endl;
 //       for(int j = 0; j< row_len; j++){
//				cout << cof[box_range(i, j)*row_len + j]<<" ";
  //      }
   //     cout<< ipiv[i] << endl;
    }
//	cout << "data2filter"<<endl;
	fftw_execute(plan_data2filter);	
	//	for(int i = 0; i < 32*32*16; i++) cout << filter[i][0]<<" "<<filter[i][1]<<"\t";
//	cout << "fft done"<<endl;
	
}


int
Filter2d::apply(double *indata, double *outdata, int len)
{
  assert(len <= size);
  int i = 0;
  for(i = 0; i < len; i++){
     data[i] = indata[i];
  }
  while(i < size) { data[i++] = 0; }
  fftw_execute(plan_data2datat);

  for(int i =0; i < size2; i++) {
	datat_filtered[i][0] = (filter[i][0]*datat[i][0] - filter[i][1]*datat[i][1]) /size;
	datat_filtered[i][1] = (filter[i][0]*datat[i][1] + filter[i][1]*datat[i][0]) /size;
  }
  fftw_execute(plan_filtered2data); /* repeat as needed */

  for(i = 0; i < len; i++){
     outdata[i] = data[i];
  }
}
