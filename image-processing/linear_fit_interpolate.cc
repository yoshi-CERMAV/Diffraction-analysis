
#include <Accelerate/Accelerate.h>
#include <valarray>
#include <iostream>
#include <fstream>
#include "Savitzky_Golay_2d.h"

using namespace std;
float *A;
int num_para = 3;
int height = 619;
int width = 487;
int length = width*height;

int fft_xn= 1024;
int fft_yn= 1024;
int data_xn = 487;
int data_yn = 619;
int shifty = (fft_yn - data_yn)/2;
int shiftx = (fft_xn - data_xn)/2;
int shiftx1 = shiftx + data_xn;

int box_x = 20;
int mid_x = box_x/2;
int box_y = 20;
int mid_y = box_y/2;
int box_size = box_x*box_y;
int skip_line = 19;

// approximate as  c0 + c1*x + c2 * y + c3 x**2 + c4 x * y + c5 **2;
int copy_to_fft(int *fdata, double *data, double *y);
int copy_back(double *datao, int *out_data);
//195 - 212, 407 - 424

int init_jacobian(float *A, int skip_line)
{
    int length = box_x * box_y;
    int mid_y = box_y/2;
    for(int i = 0; i < length; i++) A[i] = 1;
    for(int i = 0; i < box_y; i++){
        float *ptr = A + length + box_x*i;
        float *ptr1 = ptr+length;
        float *ptr2 = ptr1+length;
        float *ptr3 = ptr2+length;
        float *ptr4 = ptr3+length;
        int k = (i < mid_y )? i: i+ skip_line;
        for(int j = 0; j < box_x; j++){
            ptr[j] = j;
            ptr1[j] = k;
            ptr2[j] = j*j;
            ptr3[j] = k*j;
            ptr4[j] = k*k;
        }
    }
}

int init_jacobian(float *A)
{
    int length = box_x * box_y;
    for(int i = 0; i < length; i++) A[i] = 1;
    for(int i = 0; i < box_y; i++){
        float *ptr = A + length + box_x*i;
        float *ptr1 = ptr+length;
        float *ptr2 = ptr1+length;
        float *ptr3 = ptr2+length;
        float *ptr4 = ptr3+length;
        for(int j = 0; j < box_x; j++){
            ptr[j] = j;
            ptr1[j] = i;
            ptr2[j] = j*j;
            ptr3[j] = i*j;
            ptr4[j] = i*i;
        }
    }
}

int switch_off(float *A, int count){
    for(int i = 0; i < 6; i++) A[count+i*box_size] = 0;
}
float evaluate(float *B, int i, int j)
{
    return B[0]+ B[1]*i + B[2]*j + B[3]* i* i + B[4]*j*i + B[5]*j*j;
}


int fit_extend(double *data, int i0, int i1, int j0, int j1,int len)
{
    float *A = new float[box_size*6];
    float *B = new float[box_size];
    init_jacobian (A);
    int count = 0;
    for(int i = 0; i < box_y; i++){
        double *iptr = data + i*len;
        for(int j = 0; j < box_x; j++, count++){
                B[count] = iptr[j];
        }
    }
    int b0 = box_size;
    int b1 = 1;
    int param = 6;
    int lwork = -1;
    int info;
    float worksize;
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, &worksize, &lwork, &info);
    lwork = worksize+1;
    float *work = new float[lwork];
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, work, &lwork, &info );
    for(int j = j0; j < j1; j++){
        for(int i = i0; i < i1; i++){
            *(data+i+j*len) = evaluate(B, i, j);
        }
    }
}


int fit(int *data)
{
    float *A = new float[box_size*6];
    float *B = new float[box_size];
    init_jacobian (A);
    int count = 0;
    for(int i = 0; i < box_y; i++){
        int *iptr = data + i*width;
        for(int j = 0; j < box_x; j++, count++){
            if(iptr[j] == -1) {switch_off(A, count); B[count]= 0;}
            else{
                B[count] = iptr[j];
            }
        }
    }
    int b0 = box_size;
    int b1 = 1;
    int param = 6;
    int lwork = -1;
    int info;
    float worksize;
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, &worksize, &lwork, &info);
    lwork = worksize+1;
    float *work = new float[lwork];
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, work, &lwork, &info );
    for(int i = 0; i < box_y; i++){
        int *iptr = data + i*width;
        for(int j = 0; j < box_x; j++, count++){
            if(iptr[j] == -1) {iptr[j]=evaluate(B, j, i);}
        }
    }
}


int fill_gap(int *data, int l0)
{

    int x_range = width - box_x;
    float *A = new float[box_size*6];
    float *B = new float[box_size * x_range];
    int l1 = l0 - box_y/2;
    init_jacobian(A, skip_line);
    int *ptr = data + l1 * width;
    float *bptr = B;
    for(int col = 0; col < x_range; col++, ptr++, bptr+=box_size){
        for(int i = 0; i < box_y; i++){
            int k = (i < box_y/2)? i : i + skip_line;
            float *bptr1 = bptr + box_x *i;
            int *ptr1 = ptr + width * k;
            for(int j = 0; j < box_x; j++){
                bptr1[j] = ptr1[j];
            }
        }
    }
    int b0 = box_size;
    int b1 = x_range;
    int param = 6;
    int lwork = -1;
    int info;
    float worksize;
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, &worksize, &lwork, &info);
    lwork = worksize+1;
    float *work = new float[lwork];
    sgels_("No tranpose", &b0, &param, &b1, A, &b0, B, &b0, work, &lwork, &info );
    ptr = data+l1*width;
    
    float *C = new float[6*skip_line];
    int i1 = mid_y + skip_line;
    for(int col = 0; col < box_x; col++){
        int col1 = (col < mid_x) ? col: col+x_range;
        float *bptr = (col < mid_x) ? B : B+(x_range-1) * box_size;
        float x2 = col * col;
        for(int i = mid_y; i < i1; i++)
            ptr[i*width + col1] = bptr[0] + bptr[1] * col + bptr[2] * i + bptr[3] * x2 + bptr[4] * col * i + bptr[5] * i*i;
    }
    int mid_x2 = mid_x * mid_x;
    for(int col = 0; col < x_range; col++){
        float *bptr = B + col * box_size;
        float x2 = col*col;
        for(int i = mid_y; i < i1; i++)
            ptr[i*width + col+mid_x] = bptr[0] + bptr[1] * mid_x + bptr[2] * i + bptr[3] * mid_x2 + bptr[4] * mid_x * i + bptr[5] * i * i;
    }
    delete[] C;
    delete[] work;
    delete[] B;
    delete[] A;
}

static int *datai;
static char *buffer;

int read_file(char *filename, int *data);
int allocate_read_cbf(int w, int h);

int add_line(int *data, float &sum, float &sum2, int len, int &count)
{
    for(int i = 0; i < len; i++){
        if(data[i] == -1) continue;
        sum+= data[i];
        sum2 += data[i]*data[i];
        count++;
    }
}
int clean(int *datai, int offset_x, int offset_y, int box)
{
    int count = 0;
    int *ptr = datai+ offset_x + offset_y * width;
    int len = box*box;
    float sum = 0;
    float sum2 = 0;
    for(int i = 0; i < box; i++, ptr+=width)
        add_line(ptr, sum, sum2, box, count);
    sum /= count;
    sum2 /= count;
    float stdev = sqrt(sum2-sum);
}

int main(int argc, char *argv[])
{
    datai = new int[length];
    allocate_read_cbf(width, height);
    read_file(argv[1], datai);  //input the cbf file
    cout <<"file read"<<endl;
#ifdef AAA
    ifstream fi("hotspot");
    int i, j;
    while(fi >> i >> j ) datai[i+width*j] = -1;
    fit(datai+15+width*77);
    fit(datai+width-21+width* 64);
    fit(datai+width-21+width* 83);
    fill_gap(datai, 194);
    fill_gap(datai, 406);
    double *data = new double[1024*1024];
    double *datao = new double[1024*1024];
    cout <<"shiftx = "<<shiftx;
    for(int j = 0; j < height; j++){
        double *ptr = data + (j + shifty)*fft_xn +shiftx;
        int *ptri = datai + data_xn * j;
        for(int i = 0; i < data_xn; i++){
            ptr[i] = ptri[i];
        }
    }
    
    double *ptr =data+ shiftx + shifty * fft_xn;
    fit_extend(ptr, -shiftx, box_x, -shifty, 0, 1024);
    fit_extend(ptr + width - box_x, box_x, fft_xn - width - box_x-shiftx, -shifty, box_y, 1024);
//    fit_extend(ptr+(height-box_y)*1024, -box_x, box_x, 0, 2*box_y, 1024);
//    fit_extend(ptr + (height-box_y)*1024 +width-box_x-1, 0, box_x*2, 0, box_y*2, 1024);
    int imax = width/box_x;
    int jmax = height/box_y;
    for(int i = 1; i < imax; i++){
        fit_extend(ptr + i*box_x, 0, box_x, -shifty, 0, 1024);
//        fit_extend(ptr +(height-box_y-1)*1024 + i*box_x, 0, box_x, box_y, box_y*2, 1024);
    }
    ptr = data + shifty * fft_xn;
    double *ptr1 = ptr +width;
    int i1 = fft_xn-shiftx-width;
    int i0 = shiftx+width;
    for(int j = 0; j < height; j++){
        for(int i = 0; i < shiftx; i++){
            ptr[i + j* fft_xn] = datai[shiftx-i-1 + j*width];
        }
        for(int i = 0; i < i1; i++){
            ptr[i +i0 + j* fft_xn] = datai[width-1-i + j*width];
        }
    }
    ptr = data + (shifty+height-2) * fft_xn;
    int j1 =fft_yn - shifty-height+2;
    for(int j = 0;j < j1; j++){
        for(int i = 0; i < fft_xn; i++)
            ptr[j*fft_xn + i] = *(ptr-(j-1)*fft_xn+i);
    }
//    for(int j = 0;j < shifty ; j++){
//        for(int i = 0; i < fft_xn; i++)
//            data[j*fft_xn + i] = data[(2*shifty-j)*fft_xn + i];
//    }
//    for(int j = 1; j< jmax; j++){
//        fit_extend(ptr +j*box_y*1024, -box_x, box_x, 0, box_y, 1024);
//        fit_extend(ptr + j*box_y*1024 + width - box_x-3, 0, box_x*2, -box_y, box_y, 1024);
//    }
//#endf AAA

//    ofstream fo1("temp1.dat");
    
//    fo1.write(reinterpret_cast<char *>(data), sizeof(double)*1024*1024);
    Filter2d SGFilter (9, 9, 0, 2, 1024, 1024);
    SGFilter.apply(data, datao, 1024*1024);
//    ofstream fo2("temp2.dat");
//    fo2.write(reinterpret_cast<char *>(datao), sizeof(double)*1024*1024);
    for(int j = 0; j< height; j++){
        for(int i = 0; i < width; i++){
            datai[i + j* width] = datao[i + shiftx + (j+shifty) * 1024];
        }
    }
#endif
    ofstream fo(argv[2]);
    fo.write(reinterpret_cast<char *> (datai), length * sizeof(int));
}


int copy_to_fft(int *fdata, double *data, double *y)
{
    for(int i = 0; i < data_yn; i++){
        int *ptr = fdata+ data_xn*i;
        double *ptr1 = data + fft_xn*(i+shifty);
        for(int j = 0; j < data_xn; j++){
            ptr1[shiftx + j] = ptr[j];
        }

        for(int j = 0; j < shiftx; j++){
            double x  = j - shiftx;
            double x1 = j + data_xn;
            ptr1[j] = y[0] + x * y[1] + i * y[2] + x * i * y[3];
            ptr1[j+shiftx1] = y[0] + x1 * y[1] + i * y[2] + x1 * i * y[3];
        }
    }
    double *ptr = data;
    double *ptr1 = data + fft_xn*(shifty + data_yn);
    for(int j = 0; j < shifty; j++){
        for(int i = 0; i < fft_xn; i++, ptr++, ptr1++){
            double x = i - shiftx;
            double y0 = j - shifty;
            double y1 = j + data_yn;
            *ptr  = y[0] + x* y[1] + y0 * y[2] + x * y0 * y[3];
            *ptr1 = y[0] + x* y[1] + y1 * y[2] + x * y1 * y[3];
        }
    }
}

int copy_back(double *datao, int *out_data)
{
    for(int i = 0; i < data_yn; i++){
        int *ptr = out_data + data_xn * i;
        double *ptr1 = datao + fft_xn * (i + shifty) + shiftx;
        for(int j = 0; j < data_xn; j++){
            ptr[j] = ptr1[j];
        }
    }
}
