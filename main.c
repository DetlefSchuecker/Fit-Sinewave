/*******************************************************/
/* fit a sinewave

Detlef.Schuecker@gmx.de 03/2021

The difference equation for a sampled undamped system is

sig(n+1) = 2*cos(w)*sig(n)-sig(n-1)

This is a linear system of equations in the samples for the unknown 2*cos(w).
I solve that overdetermined system using linear regression.

I combine samples which are close to half the nyquist rate, i.e. 4 samples per sinewave.
cos(w) is flat at the top for small w, so we run in precision problems for
small w and we thus have to shuffle the samples, see variable ddoff in the code.

*/
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/*******************************************************/
uint32_t fitsinus(int16_t  * sig,  // the signal
                  int32_t n,      // length
                  double * amp,   // returned amplitude
                  double * wn,   //  returned frequency
                  double * ph,  //  returned phase
                  double * offs,//  returned offset
                  double * err  //  diff. reconstr./original in dB
                  )
/*******************************************************/
{
   double rere,imim,offf;
   double a11,a12,a13;
   double a21,a22,a23;
   double a31,a32,a33;
   double det,sa,sb;
   double coffa,coffb;
   double s1,s2,oc;
   double aa,bb,cc,dd,ee;
   double ff,gg,hh,ii;

   int32_t k,istate,icnt;
   int32_t iedge,inew,ddoff;

   // calc. coarse offset
   oc =0.0;
   for(k=0;k<n;k++) oc += sig[k];
   oc = oc/n;

   //printf("oc %.2f \n", oc);

   // do coarse frequency, check for zero crossings
   // check for bouncing close to zero
   istate=((sig[0]-oc)<0.0);
   icnt  =0;
   iedge =0; // number of zero crossings

   for(k=1;k<n;k++){
     inew=((sig[k]-oc)<0.0);
     if(istate==1) {
       if(inew==0){icnt++; if(icnt>7) {iedge++;istate = 0;}} else icnt=0;
     }else{
       if(inew==1){icnt++; if(icnt>7) {iedge++;istate = 1;}} else icnt=0;
     }
   }

   if(iedge==0) return -1; // no zero crossing ?
   ddoff=n/(2*iedge); // distance to get near nyquist/2

   //printf("iedge %d \n",iedge);
   //printf("ddoff %d \n",ddoff);

   // define two signals
   // s1=sig[0..n-1-2*ddoff]+sig[2*ddoff..n-1], length n-2*ddoff
   // s2=sig[ddoff..n-1-ddoff]                , length n-2*ddoff
   // and a matrix
   // M =[s2 1]

   // calc. M'*M =[[aa bb];[bb cc]]
   // and   M'*s2 = [dd ee]'

   aa=0.0;
   bb=0.0;
   cc=n-2*ddoff;
   dd=0.0;
   ee=0.0;

   for(k=0;k<n-2*ddoff-1;k++){
      s1=sig[k   ]+sig[k+2*ddoff];
      s2=sig[k+ddoff];
      aa += s2*s2;
      bb += s2;
      dd += s2*s1;
      ee += s1;
   }

   det=aa*cc-bb*bb;
   if(det==0.0) return -2; // no solution ?
   det=1.0/det;

   coffa = det*( cc*dd-bb*ee);
   coffb = det*(-bb*dd+aa*ee);

   wn[0]=acos(coffa/2.0)/ddoff;

   //printf("wnr %.8f\n",wn[0]);

   // now mix with quadrature carrier plus offset
   // matrix M   =[sa sb 1] = [cos(wn) sin(wn) 1], length n
   // matrix M'*M=[[aa bb cc];[bb dd ee];[cc ee ff]], length 3x3

   aa=0.0;bb=0.0;cc=0.0;
   dd=0.0;ee=0.0;ff=0.0;
   gg=0.0;hh=0.0;ii=0.0;

   for(k=0;k<n-1;k++){
      sa=cos(k*wn[0]);
      sb=sin(k*wn[0]);

      aa+= sa*sa;
      bb+= sa*sb;
      cc+= sa;
      dd+= sb*sb;
      ee+= sb;
      ff+= 1.0; // :)

      gg+= sa*sig[k];
      hh+= sb*sig[k];
      ii+=    sig[k];
   }

   det=aa*(dd*ff-ee*ee)+bb*(cc*ee-bb*ff)+cc*(bb*ee - cc*dd);
   if(det==0.0) return -3; // no solution ?
   det=1.0/det;

   a11=dd*ff-ee*ee;
   a21=cc*ee-bb*ff;
   a31=bb*ee-cc*dd;

   a12=cc*ee-bb*ff;
   a22=aa*ff-cc*cc;
   a32=bb*cc-aa*ee;

   a13=bb*ee-cc*dd;
   a23=bb*cc-aa*ee;
   a33=aa*dd-bb*bb;

   //printf("11 %.8f %.8f %.8f\n",a11,a12,a13);
   //printf("22 %.8f %.8f %.8f\n",a21,a22,a23);
   //printf("33 %.8f %.8f %.8f\n",a31,a32,a33);

   rere=det*(a11*gg+a12*hh+a13*ii);
   imim=det*(a21*gg+a22*hh+a23*ii);
   offf=det*(a31*gg+a32*hh+a33*ii);

   //printf("rere %.8f\n",rere);
   //printf("imim %.8f\n",imim);
   //printf("offf %.8f\n",offf);

   amp[0] = sqrt(rere*rere+imim*imim);
   ph[0]  = atan2(-imim,rere);
   offs[0]= offf;

   //printf("amp %.8f\n",amp[0]);
   //printf("ph %.8f\n",ph[0]);
   //printf("offs %.8f\n",offs[0]);

   // reconstruct signal
   cc=0.0;
   ff=0.0;
   for(k=0;k<n;k++){
      aa=amp[0]*cos(k*wn[0]+ph[0]);
      ee=(sig[k]-offs[0]);
      ff+=ee*ee;
      dd=aa-ee;
      cc+=dd*dd;
      sig[k]=round(aa+offs[0]); // this line overwrites original
                               // comment out if necessary
   }

   err[0]=10.0f*log10(ff/(cc+1E-50));

   //printf("err %.2f\n",err[0]);

   return 0;
}

/*******************************************************/
int main()
/*******************************************************/
{ double amp ,wn ,ph ,offs ;
  double ampr,wnr,phr,offsr,err;

  int16_t sig[100000];
  int32_t n,k,ret;


n=6000;
wn= 2*3.14*100.1/n;
amp=8345.0;
ph=0.56789;
offs=0.0;

for(k=0;k<n;k++)
 sig[k]=(int16_t)(amp*cos(k*wn+ph)+offs+0.0*((rand()%201)-100));

//for(k=0;k<n;k++)
// printf("%d %d\n",k,sig[k]);

ret=fitsinus(sig,n,&ampr,&wnr,&phr,&offsr,&err);

printf("amp  %.8f\r\n",amp-ampr);
printf("wn  %.8f\r\n",wn-wnr);
printf("ph  %.8f\r\n",ph-phr);
printf("off  %.8f\r\n",offs-offsr);
printf("err  %.8f\r\n",err);

printf("Hello world!\n\r %d \n\r",
       ret);

return 0;
}
