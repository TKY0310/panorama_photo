#include"image.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define MAX 30
#define TRIAL 1000

double VP(double*a,double*b,int N){
  double s=0;
  int i;
  for(i=0;i<N;i++) s += a[i] * b[i] ;
  return s;
}

void VSS(double*d,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] *= s;
}

void VSA(double*d,double*a,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] += a[i] * s;
}

void MatrixCopy(Matrix*mtD,Matrix*mt){
  memmove(mtD->data,mt->data,mt->W*mt->H*sizeof(double));
}

void MatrixCopyT(Matrix*mtD,Matrix*mt){
  int i,j;
  for(i=0;i<mtD->H;i++)
    for(j=0;j<mtD->W;j++)
      Elem(mtD,i,j) = Elem(mt,j,i);
}

void MatrixMultT(Matrix*mtD,Matrix*mtA,Matrix*mtB){
  // D = A B^T
  int i,j;
  for(i=0;i<mtA->H;i++)
    for(j=0;j<mtB->H;j++)
      Elem(mtD,i,j) = VP( Row(mtA,i), Row(mtB,j), mtA->W);
}

void MatrixQRDecompColMajor(Matrix*mtR,Matrix*mt){
  // Gram-Schmidt orthonormalization (R and Q)
  double t, *aT[] = { Row(mt,0), Row(mt,1), Row(mt,2), Row(mt,3),Row(mt,4),Row(mt,5),Row(mt,6),Row(mt,7)} ;
  int W = mt->W;
  MatrixClear(mtR);
	int i,j;
	for(i = 0; i<mt->H;i++){
		for(j = 0;j!=i;j++){
			Elem(mtR,j,i) = t = VP(aT[j], aT[i], W);
  			VSA(aT[i], aT[j], -t, W);
		}
		Elem(mtR,i,i) = t = sqrt(VP(aT[i],aT[i],W));
  		VSS(aT[i], 1/t, W);
	}
}

void MatrixSimeqLr(Matrix*mtB,Matrix*mtR){
  // B = B L^{-1}
  double * B = Row(mtB,0);
	int i,j;
	for(i = mtB->W-1;i>=0;i--){
		for(j = i+1; j<mtB->W;j++){
			B[i] -= B[j]*Elem(mtR,i,j);
		}
		B[i] /= Elem(mtR,i,i);
	}
}

void ImageFeature(Matrix*im2,Image*im){
  int x,y,u,v,W=7,ix,iy;
  for(y=W+1;y<im->H-W-1;y++) for(x=W+1;x<im->W-W-1;x++){
    double ixx,ixy,iyy;
    ixx=iyy=ixy=0;
    for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
      ix=IElem(im, x+u+1, y+v, 1) - IElem(im, x+u-1, y+v, 1);
      iy=IElem(im, x+u,  y+v+1,1) - IElem(im, x+u,  y+v-1,1);
    	
      ixx+=ix*ix; // 
      iyy+=iy*iy;
      ixy+=ix*iy;
    }
  	double b = -(ixx+iyy),c = ixx*iyy-(ixy*ixy);
  	double ramda = (-b - sqrt(b*b - 4*c)) /2;
  	
    DElem(im2,x,y) = ramda; //  
  }
}

int MatrixLocalMax(int w[][2], double k[], Matrix*im2){
  int x,y,u,v,W=7,n=0,a;
  for(y=W+1;y<im2->H-W-1;y++) for(x=W+1;x<im2->W-W-1;x++){
    double max=-1;
    for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
    	if(DElem(im2,x+u,y+v)>max) max =  DElem(im2,x+u,y+v);
    }
  	if(max == DElem(im2,x,y)){
  		a = n;
  		if(n < MAX) n++;
  		for(;a>0 && k[a-1] < DElem(im2,x,y);a--){
  			w[a][0] = w[a-1][0];
  			w[a][1] = w[a-1][1];
  			k[a] = k[a-1];
  		}
  			w[a][0] = x;
  			w[a][1] = y;
  			k[a] = DElem(im2,x,y);
  	}
  }
  return n; //
}

double ImageSSD(Image*im,int x1,int y1, Image*im2,int x2,int y2){
  int i,j,W=7;
  double sr=0,sg=0,sb=0,dr,dg,db;
  for(i=-W;i<=W;i++) for(j=-W;j<=W;j++){
    dr  = IElem(im, x1+j, y1+i, 0) - IElem(im2, x2+j , y2+i, 0);
    dg  = IElem(im, x1+j, y1+i, 1) - IElem(im2, x2+j , y2+i, 1);
    db  = IElem(im, x1+j, y1+i, 2) - IElem(im2, x2+j , y2+i, 2);
    sr += dr*dr;
    sg += dg*dg;
    sb += db*db;
  }
  return sr+sg+sb;
}


void calcSSDtable(Matrix*mt,
		  Image*im ,int x1[][2],int N1,
		  Image*im2,int x2[][2],int N2){
  int i,j;
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      Elem(mt,i,j) = ImageSSD(im ,x1[i][0],x1[i][1],
			      im2,x2[j][0],x2[j][1]);
}

int matchMethod2(double w[][4],Matrix*mt,
                 Image*im ,int x1[][2],int N1,
                 Image*im2,int x2[][2],int N2){
  int i,j,k,vali,valj,n=0;
   for(n=0;n<N1;n++){
   	double sm=INFINITY;
  	for(i=0;i<N1;i++)for(j=0;j<N2;j++){
      double t;
  	  t=Elem(mt,i,j);
      if(sm>t) sm=t, vali=i,valj=j;
    }
    w[n][0] = x1[vali][0];
   	w[n][1] = x1[vali][1];
   	w[n][2] = x2[valj][0];
   	w[n][3] = x2[valj][1];
   	
    for(k=0;k<N1;k++) Elem(mt,k,valj)=INFINITY;
    for(k=0;k<N2;k++) Elem(mt,vali,k)=INFINITY;
  }
 
  return n;
}

void calcHomography(double H[][3],double w[][4],int rndAry[], Matrix*cmA,Matrix*vt,Matrix*mtR,Matrix*tmp){
	int i;
	int a=rndAry[0], b=rndAry[1], c=rndAry[2], d=rndAry[3];
	
    double ww[][4]={
    w[a][0], w[a][1], w[a][2], w[a][3],
    w[b][0], w[b][1], w[b][2], w[b][3],    
    w[c][0], w[c][1], w[c][2], w[c][3],
  	w[d][0], w[d][1], w[d][2], w[d][3],
  };
	for(i=0;i<4;i++){
    Elem(cmA,0,i*2  )=ww[i][0];
    Elem(cmA,1,i*2  )=ww[i][1];
    Elem(cmA,2,i*2  )=1;
    Elem(cmA,3,i*2  )=0;
    Elem(cmA,4,i*2  )=0;
    Elem(cmA,5,i*2  )=0;
    Elem(cmA,6,i*2  )=-ww[i][0]*ww[i][2];
    Elem(cmA,7,i*2  )=-ww[i][1]*ww[i][2];
    Elem(cmA,0,i*2+1)=0;
    Elem(cmA,1,i*2+1)=0;
    Elem(cmA,2,i*2+1)=0;
    Elem(cmA,3,i*2+1)=ww[i][0];
    Elem(cmA,4,i*2+1)=ww[i][1];
    Elem(cmA,5,i*2+1)=1;
    Elem(cmA,6,i*2+1)=-ww[i][0]*ww[i][3];
    Elem(cmA,7,i*2+1)=-ww[i][1]*ww[i][3];
    Elem(vt ,0,i*2  )=ww[i][2];
    Elem(vt ,0,i*2+1)=ww[i][3];
  }
	
	MatrixQRDecompColMajor(mtR,cmA);
	MatrixMultT(tmp,vt,cmA);
	MatrixSimeqLr(tmp,mtR);
	
	H[0][0] = Elem(tmp,0,0);
	H[0][1] = Elem(tmp,0,1);
	H[0][2] = Elem(tmp,0,2);
	H[1][0] = Elem(tmp,0,3);
	H[1][1] = Elem(tmp,0,4);
	H[1][2] = Elem(tmp,0,5);
	H[2][0] = Elem(tmp,0,6);
	H[2][1] = Elem(tmp,0,7);
	H[2][2] = 1;
}

void initRndAry(int rndAry[]){
	int i;
    for(i=0;i<MAX;i++) rndAry[i]=i;
    srand(__rdtsc());
}

void chooseFourNumbers(int rndAry[]){
	int i;
  for(i=0;i<4;i++){
  	int j, t;
  	j=(int)((long long)rand()*(MAX-i)/(RAND_MAX+1LL))+i; 
  	t=rndAry[i]; rndAry[i]=rndAry[j]; rndAry[j]=t;
  }
}

int calcScore(double H[][3],double w[][4]){
	int i,n=0;
	
	for(i=0;i<MAX;i++){
    double x=w[i][0], y=w[i][1],
	       u=w[i][2], v=w[i][3],
	x_prime = H[0][0]*x + H[0][1]*y + H[0][2],
	y_prime = H[1][0]*x + H[1][1]*y + H[1][2], 
	z_prime = H[2][0]*x + H[2][1]*y + H[2][2], 
	du = u - (x_prime / z_prime),
	dv = v - (y_prime / z_prime); 
		if(du*du+dv*dv<5) n++;
  }
	return n;
}

void Match(double w[][4],Image*im,Image*im2,int x1[][2],int x2[][2],int N1,int N2){
	int nm;
	Matrix*mt;
	mt = MatrixAlloc(N1,N2);
 
	calcSSDtable(mt,im,x1,N1,im2,x2,N2);
	nm = matchMethod2(w,mt,im,x1,N1,im2,x2,N2);
}

void ransac(double w[][4],double bestH[][3]){
	int rndAry[MAX];
	int trial,score,best_score =0;
	Matrix*cmA,*vt,*mtR,*tmp; 
	cmA=MatrixAlloc(8,8);
	vt=MatrixAlloc(1,8);
	mtR=MatrixAlloc(8,8);
	tmp=MatrixAlloc(1,8);
	
	initRndAry(rndAry);
	
	 for(int trial=0;trial<TRIAL;trial++){
   double H[3][3]; 
   chooseFourNumbers(rndAry);
   calcHomography(H,w,rndAry, cmA,vt,mtR,tmp);
   score=calcScore(H,w);
   if(score > best_score){ 
   	for(int i=0; i<3; i++){
   		for(int j=0; j<3; j++){
   			bestH[i][j] = H[i][j];
   			best_score = score;
   		}
	}
   }
 }
}

int FeatureExtraction(Image*im,int x[][2],double valx[]){
	int n;
	Matrix*imf;
	imf = MatrixAlloc(im->H,im->W);
	ImageFeature(imf,im);
	n = MatrixLocalMax(x,valx,imf);
	return n;
}

void createPanorama(Image*id,Image*im,Image*im2,double bestH[][3],double alpha){
	double m0d[][3]={
      1,0,-100,
      0,1,-100,
      0,0,1
    },m1d[3][3];
	mult33(m1d,bestH,m0d);
	ImageImageProjectionAlpha(id,im2,m1d,alpha);
}

int main(int ac,char**av){
	double bestH[3][3];
	double m0d[][3]={
      1,0,-100,
      0,1,-100,
      0,0,1
    };
	Image *im,*im2,*imout;
	double alpha;
	
	int x1[MAX+1][2],x2[MAX+1][2],N1,N2;
	double w[MAX+1][4],valx1[MAX+1],valx2[MAX+1];

	imout=ImageAlloc(1024,768);
	ImageClear(imout);
	
	im = ImageRead(av[1]);
	N1 = FeatureExtraction(im,x1,valx1);
	
	im2 = ImageRead(av[2]);
	N2 = FeatureExtraction(im2,x2,valx2);

	Match(w,im,im2,x1,x2,N1,N2);

	ransac(w,bestH);
	
	alpha = 1 / (double)(ac-1);
	ImageImageProjectionAlpha(imout,im,m0d,alpha);
	createPanorama(imout,im,im2,bestH,alpha);
	
	if(ac>3){
	Image *im3;
	int x3[MAX+1][2],N3;
	double valx3[MAX+1];
	im3 = ImageRead(av[3]);
	N3 = FeatureExtraction(im3,x3,valx3);
	
	Match(w,im,im3,x1,x3,N1,N3);

	ransac(w,bestH);
	createPanorama(imout,im,im3,bestH,alpha);
	
	}
	 ImageWrite("test2.jpg",imout);
	return 0;
}