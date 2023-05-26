#include "complex_math.h"

gsl_complex gsl_complex_div_real(gsl_complex a, double b)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a)/b,GSL_IMAG(a)/b);
    
    return z;
}


gsl_complex gsl_complex_conjugate(gsl_complex a)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a),-GSL_IMAG(a));
    
    return z;
}

gsl_complex gsl_complex_mul(gsl_complex a, gsl_complex b)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a)*GSL_REAL(b)-GSL_IMAG(a)*GSL_IMAG(b),GSL_REAL(a)*GSL_IMAG(b)+GSL_IMAG(a)*GSL_REAL(b));
    
    return z;
}

gsl_complex gsl_complex_add(gsl_complex a, gsl_complex b)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a)+GSL_REAL(b),GSL_IMAG(a)+GSL_IMAG(b));
    
    return z;
}

gsl_complex gsl_complex_sub(gsl_complex a, gsl_complex b)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a)-GSL_REAL(b),GSL_IMAG(a)-GSL_IMAG(b));
    
    return z;
}

gsl_complex gsl_complex_negative(gsl_complex a)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,-GSL_REAL(a),-GSL_IMAG(a));
    
    return z;
}

gsl_complex gsl_complex_sub_real (gsl_complex a, double x)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z,GSL_REAL(a)-x,GSL_IMAG(a));
    
    return z;
}

gsl_complex gsl_complex_rect (double x, double y)
{
    gsl_complex z;
    GSL_SET_COMPLEX(&z, x, y);
    
    return z;
}

gsl_complex gsl_complex_inverse (gsl_complex x)
{
    gsl_complex z;
    double r;
    r = gsl_complex_abs2(x);
    GSL_SET_COMPLEX(&z, GSL_REAL(x)/r, -GSL_IMAG(x)/r);
    
    return z;
}

double gsl_complex_abs2(gsl_complex a)
{
    return GSL_REAL(a)*GSL_REAL(a)+GSL_IMAG(a)*GSL_IMAG(a);
}

double gsl_complex_abs(gsl_complex a)
{
    return sqrt(gsl_complex_abs2(a));
}

int gsl_complex_eq(gsl_complex a, gsl_complex b)
{
    if( GSL_REAL(a)==GSL_REAL(b)&& GSL_IMAG(a)==GSL_IMAG(b))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

size_t  gsl_vector_complex_max_index(const gsl_vector_complex *v)
{
    size_t idx;
    gsl_complex max, element;
    double absMax;
 
    idx = 0;
    max = gsl_vector_complex_get(v, idx);
    absMax = gsl_complex_abs(max);
    for (size_t i = 1; i < v->size; i++)
    {
        element = gsl_vector_complex_get(v, i);
        // Test if new maximum and update
        if (absMax < gsl_complex_abs(element))
        {
            max = element;
            absMax = gsl_complex_abs(max);
            idx = i;
        }
    } 
   return idx;
}
