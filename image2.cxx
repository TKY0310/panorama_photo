extern "C" {
#include"image.h"
}

#include<windows.h>
#include<gdiplus.h>
#pragma comment(lib,"gdiplus.lib")

class _GdiPlus {
  Gdiplus::GdiplusStartupInput gdiplusStartupInput;
  ULONG_PTR gdiplusToken;
public:
  _GdiPlus(){   Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL); }
  ~_GdiPlus(){  Gdiplus::GdiplusShutdown(gdiplusToken); }
} _gdiplus;

#define upto4(a) (((a)+3)&-4)

Image*ImageRead(const char*name){
  wchar_t wname[MAX_PATH];
  int len=MultiByteToWideChar( 932, MB_PRECOMPOSED, name,-1, wname,MAX_PATH );

  Gdiplus::Bitmap img(wname,1);
  int wdt, hgt, BPP=3;
  Gdiplus::Rect rc(0,0,wdt=img.GetWidth(),hgt=img.GetHeight());
  Gdiplus::BitmapData bd;
  Gdiplus::Status st=img.LockBits(&rc, Gdiplus::ImageLockModeRead, PixelFormat24bppRGB, &bd);

  Image*im=ImageAlloc(wdt,hgt);
  BYTE*wp=im->data, *rp;
  for(int y=0;y<hgt;y++){
    rp=(BYTE*)bd.Scan0+(y)*bd.Stride; // hgt-1-
    for(int x=0;x<wdt;x++,wp+=BPP,rp+=3)
      wp[0]=rp[2], wp[1]=rp[1], wp[2]=rp[0];
//    for(int y=0;y<im->H;y++) fread(im->data+im->W*3*(im->H-1-y),1,im->W*3,fp);
  }
  img.UnlockBits(&bd);
  return im;
}

void ImageWrite(const char*name,Image*im){
  int wdt=im->W, hgt=im->H, pitch=upto4(wdt*3);
  BYTE*img=new BYTE[hgt*pitch], *rp, *wp;
  int x,y;
  for(y=0;y<hgt;y++){
    rp=im->data+(y)*wdt*3; // hgt-1-
    wp=img+y*pitch;
    for(x=0;x<wdt;x++,wp+=3,rp+=3)
      wp[0]=rp[2], wp[1]=rp[1], wp[2]=rp[0];
  }
  
  const static CLSID encoder[]={
    0x557cf401, 0x1a04, 0x11d3, {0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e},
  };

  Gdiplus::Bitmap b(wdt,hgt,pitch,PixelFormat24bppRGB, img);
  WCHAR wname[MAX_PATH];
  MultiByteToWideChar(CP_ACP, 0, name, -1, wname, MAX_PATH);
  Gdiplus::Status hr=b.Save(wname, encoder, NULL);
  if(hr!=Gdiplus::Ok) fprintf(stderr,"ERROR in ImageWrite %s\n",name);
  delete[] img;
}
