
#ifndef PRECISION_H
#define PRECISION_H

//#include <cmath> 
#include <math.h>

#ifdef SINGLE_PRECISION
       typedef float RealVar;
      
       #define Fabs(x) ( fabsf(x) ) 
       #define Sin(x) ( sinf(x) ) 
       #define Cos(x) ( cosf(x) ) 
       #define Acos(x) ( acosf(x) ) 
       #define Atan(x) ( atanf(x) ) 
       #define Atan2(y,x) ( atan2f(y,x) ) 
       #define Exp(x) ( expf(x) ) 
       #define Log(x) ( logf(x) )
       #define Erfc(x) ( erfcf(x) )
       #define Sqrt(x) ( sqrtf(x) ) 
       #define Floor(x) ( floorf(x) ) 
       #define Modf(x,y) ( modff(x,y) ) 
#else
       typedef double RealVar;

       #define Fabs(x) ( fabs(x) ) 
       #define Sin(x) ( sin(x) ) 
       #define Cos(x) ( cos(x) ) 
       #define Acos(x) ( acos(x) ) 
       #define Atan(x) ( atan(x) ) 
       #define Atan2(y,x) ( atan2(y,x) ) 
       #define Exp(x) ( exp(x) )
       #define Log(x) ( log(x) )
       #define Erfc(x) ( erfc(x) )
       #define Sqrt(x) ( sqrt(x) ) 
       #define Floor(x) ( floor(x) ) 
       #define Modf(x,y) ( modf(x,y) ) 
#endif

#endif

