#include "GeodesicEquationSolver.h"



double Dnu(const double nu, const Coord& x, const Tensor3& v, Metric& metric)
{
    // Vincent, Gourgoulhon, Novak 2012
    double alpha = metric.GetAlpha(x);
    Tensor3 dl_alpha = metric.GetDerivAlpha_l(x);
    Tensor3x3 K_ll = metric.GetExtrCurv_ll(x);
    double dnu = 0;
    for(int i=1; i<4; i++)
    {
        dnu -= v[i] * dl_alpha[i];
        for(int j=1; j<4; j++)
            dnu += alpha * K_ll[{i,j}] * v[i]*v[j];
    }
    dnu *= nu;
    return dnu;
}
Tensor3 Dx(const Coord& x, const Tensor3& v, Metric& metric)
{
    // Vincent, Gourgoulhon, Novak 2012     evolves 3 velocity measured by eulerian observer v^mu, p^mu = nu(n^mu+v^mu)
    double alpha = metric.GetAlpha(x);
    Tensor3 beta_u = metric.GetBeta_u(x);
    Tensor3 dx(0.0);
    for(int k=1; k<4; k++)
        dx[k] = alpha*v[k] - beta_u[k];
    return dx;
}
Tensor3 Dv(const Coord& x, const Tensor3& v, Metric& metric)
{
    // Vincent, Gourgoulhon, Novak 2012     evolves 3 velocity measured by eulerian observer v^mu, p^mu = nu(n^mu+v^mu)
    double alpha = metric.GetAlpha(x);
    Tensor3 dl_Alpha = metric.GetDerivAlpha_l(x);
    Tensor3x3 dl_Beta_u = metric.GetDerivBeta_lu(x);
    Tensor3x3 K_ll = metric.GetExtrCurv_ll(x);
    Tensor3x3 gamma_uu = metric.GetGamma_uu(x);
    Tensor3x3x3 dl_gamma_ll = metric.GetDerivGamma_lll(x);
    Tensor3 dv(0.0);
    for(int k=1; k<4; k++)
    for(int i=1; i<4; i++)
    {
        dv[k] += v[k]*dl_Alpha[i]*v[i] - gamma_uu[{k,i}]*dl_Alpha[i] - v[i]*dl_Beta_u[{i,k}];
        for(int j=1; j<4; j++)
        {
            dv[k] += 2.0*alpha*gamma_uu[{i,k}]*K_ll[{i,j}]*v[j] - alpha*v[k]*K_ll[{i,j}]*v[i]*v[j];
            for(int l=1; l<4; l++)
                dv[k] -= alpha*0.5*gamma_uu[{k,l}] * (dl_gamma_ll[{i,l,j}] + dl_gamma_ll[{j,i,l}] - dl_gamma_ll[{l,i,j}]) * v[i]*v[j];
        }
    }
    return dv;
}

template<int direction, class Coord>
double Euler_GeodesicEquation(double dt, Coord& x, Tensor3& v, Metric& metric)
{
    // nu0 = 1, s = nu/nu0 = nu.
    int N = 2;
    double ddt = dt/N;
    double nu = 1;
    for(int i=0; i<N; i++)
    {
        double dnu = Dnu(1,x,v,metric);
        Tensor3 dx = Dx(x,v,metric);
        Tensor3 dv = Dv(x,v,metric);
        for(int k=1; k<4; k++)
        {
            x[k] += dx[k]*ddt*direction;
            v[k] += dv[k]*ddt*direction;
        }
        nu += ddt*dnu;
    }
    return nu;
}

constexpr double maxTE = 1e-6;
double RK45_GeodesicEquation(double dt, Coord& x, Tensor3& v, Metric& metric)
{
    Tensor3 k1x, k2x, k3x, k4x, k5x, k6x;
    Tensor3 k1v, k2v, k3v, k4v, k5v, k6v;
    double k1nu, k2nu, k3nu, k4nu, k5nu, k6nu;
    Coord xTemp;
    Tensor3 vTemp;
    double nuTemp;
    Tensor3 dx;
    Tensor3 dv;
    double dnu;

    // nu0 = 1, s = nu/nu0 = nu.
    double nu = 1;
    double totalTime = 0;
    double ddt = metric.grid.dt;
    while(totalTime < dt)
    {
        ddt = std::min(ddt,dt-totalTime);
        k1x = Dx(x,v,metric);
        k1v = Dv(x,v,metric);
        k1nu= Dnu(nu,x,v,metric);
        for(int d=1; d<4; d++)
        {
            xTemp[d] = x[d] + ddt*( (1.0/4.0)*k1x[d] );
            vTemp[d] = v[d] + ddt*( (1.0/4.0)*k1v[d] );
        }
        nuTemp = nu + ddt*( (1.0/4.0)*k1nu );

        k2x = Dx(xTemp,vTemp,metric);
        k2v = Dv(xTemp,vTemp,metric);
        k2nu= Dnu(nuTemp,xTemp,vTemp,metric);
        for(int d=1; d<4; d++)
        {
            xTemp[d] = x[d] + ddt*( (3.0/32.0)*k1x[d] + (9.0/32.0)*k2x[d] );
            vTemp[d] = v[d] + ddt*( (3.0/32.0)*k1v[d] + (9.0/32.0)*k2v[d] );
        }
        nuTemp = nu + ddt*( (3.0/32.0)*k1nu + (9.0/32.0)*k2nu );

        k3x = Dx(xTemp,vTemp,metric);
        k3v = Dv(xTemp,vTemp,metric);
        k3nu= Dnu(nuTemp,xTemp,vTemp,metric);
        for(int d=1; d<4; d++)
        {
            xTemp[d] = x[d] + ddt*( (1932.0/2197.0)*k1x[d] - (7200.0/2197.0)*k2x[d] + (7296.0/2197.0)*k3x[d] );
            vTemp[d] = v[d] + ddt*( (1932.0/2197.0)*k1v[d] - (7200.0/2197.0)*k2v[d] + (7296.0/2197.0)*k3v[d] );
        }
        nuTemp = nu + ddt*( (1932.0/2197.0)*k1nu - (7200.0/2197.0)*k2nu + (7296.0/2197.0)*k3nu );

        k4x = Dx(xTemp,vTemp,metric);
        k4v = Dv(xTemp,vTemp,metric);
        k4nu= Dnu(nuTemp,xTemp,vTemp,metric);
        for(int d=1; d<4; d++)
        {
            xTemp[d] = x[d] + ddt*( (439.0/216.0)*k1x[d] - (8.0)*k2x[d] + (3680.0/513.0)*k3x[d] - (845.0/4104.0)*k4x[d] );
            vTemp[d] = v[d] + ddt*( (439.0/216.0)*k1v[d] - (8.0)*k2v[d] + (3680.0/513.0)*k3v[d] - (845.0/4104.0)*k4v[d] );
        }
        nuTemp = nu + ddt*( (439.0/216.0)*k1nu - (8.0)*k2nu + (3680.0/513.0)*k3nu - (845.0/4104.0)*k4nu );

        k5x = Dx(xTemp,vTemp,metric);
        k5v = Dv(xTemp,vTemp,metric);
        k5nu= Dnu(nuTemp,xTemp,vTemp,metric);
        for(int d=1; d<4; d++)
        {
            xTemp[d] = x[d] + ddt*( -(8.0/27.0)*k1x[d] + (2.0)*k2x[d] - (3544.0/2565.0)*k3x[d] + (1859.0/4104.0)*k4x[d] - (11.0/44.0)*k5x[d] );
            vTemp[d] = v[d] + ddt*( -(8.0/27.0)*k1v[d] + (2.0)*k2v[d] - (3544.0/2565.0)*k3v[d] + (1859.0/4104.0)*k4v[d] - (11.0/44.0)*k5v[d] );
        }
        nuTemp = nu + ddt*( -(8.0/27.0)*k1nu + (2.0)*k2nu - (3544.0/2565.0)*k3nu + (1859.0/4104.0)*k4nu - (11.0/44.0)*k5nu );

        k6x = Dx(xTemp,vTemp,metric);
        k6v = Dv(xTemp,vTemp,metric);
        k6nu= Dnu(nuTemp,xTemp,vTemp,metric);
        for(int d=1; d<4; d++)
        {
            dx[d] = (16.0/135.0)*k1x[d] + (6656.0/12825.0)*k3x[d] + (28561.0/56430.0)*k4x[d] - (9.0/50.0)*k5x[d] + (2.0/55.0)*k6x[d];
            dv[d] = (16.0/135.0)*k1v[d] + (6656.0/12825.0)*k3v[d] + (28561.0/56430.0)*k4v[d] - (9.0/50.0)*k5v[d] + (2.0/55.0)*k6v[d];
        }
        dnu = (16.0/135.0)*k1nu + (6656.0/12825.0)*k3nu + (28561.0/56430.0)*k4nu - (9.0/50.0)*k5nu + (2.0/55.0)*k6nu;

        // Truncation Error (TE) estimate:
        double TE = 0;
        for(int d=1; d<4; d++)
            TE += pow( ddt*((1.0/360.0)*k1v[d] - (128.0/4275.0)*k3v[d] - (2197.0/75240.0)*k4v[d] + (1.0/50.0)*k5v[d] + (2.0/55.0)*k6v[d]) , 2);
        TE = sqrt(TE);

        if(TE > maxTE)
        {// Reject new point and redo step:
            ddt *= 0.95 * pow(maxTE / TE, 1.0/4.0);
            continue;
        }
        else
        {// Accept new point:
            for(int d=1; d<4; d++)
            {
                x[d] += ddt*dx[d];
                v[d] += ddt*dv[d];
            }
            nu += ddt*dnu;
            
            totalTime += ddt;
            ddt *= 0.95 * pow(maxTE / (TE + 1e-50), 1.0/4.0);

            if(metric.InsideBH(x) or metric.grid.OutsideDomain(x))
                break;
        }
    }
    return nu;
}