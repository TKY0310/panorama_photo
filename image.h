#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct _Image {
  unsigned char*data;
  int W,H;
} Image;

typedef struct _Matrix{
  double *data;
  int W,H;
} Matrix;

#define IElem(_im,_x,_y,_c) (_im) -> data[(_y)*(_im)->W*3 +(_x)*3+(_c)]
#define isInsideImage(is,u,v) (0<=u && u<is->W && 0<=v && v<is->H)
#define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
#define DElem(_a,_b,_c) (_a)->data[(_a)->W*(_c)+(_b)]
#define Row(_a,_b)     ((_a)->data+(_a)->W*(_b))

Image*ImageAlloc(int W,int H);
Image*ImageRead(const char*name);
void ImageFree(Image*im);
void ImageWrite(const char*name,Image*im);
void ImageClear(Image*im);
void ImageImageProjectionAlpha(Image*id,Image*is,double a[3][3],double alpha);
void ImageDrawBox(Image*im,int x,int y);
void mult33(double d[3][3],double a[3][3],double b[3][3]);
Matrix*MatrixAlloc(int _H,int _W);
void MatrixClear(Matrix*mt);