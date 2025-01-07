#include "likelihood.h"

void Swap(double* a, double* b) 
{ 
  double t = *a; 
  *a = *b; 
  *b = t; 
  return;
}

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
double Partition (double arr[], int low, int high) 
{ 
  double pivot = arr[high]; // pivot 
  int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far

  for (int j = low; j <= high - 1; j++) 
  { 
    // If current element is smaller than the pivot 
    if (arr[j] < pivot) 
    { 
      i++; // increment index of smaller element 
      Swap(&arr[i], &arr[j]); 
    } 
  } 
  Swap(&arr[i + 1], &arr[high]); 
  return (i + 1); 
}

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void QuickSort(double arr[], int low, int high) 
{ 
  if (low < high) 
  { 
    /* pi is partitioning index, arr[p] is now 
    at right place */
    int pi = Partition(arr, low, high); 

    // Separately sort elements before 
    // partition and after partition 
    QuickSort(arr, low, pi - 1); 
    QuickSort(arr, pi + 1, high); 
  }
}

/*
    CAUTION: Remove median must be called on each
    sector data separately
*/
void Remove_Median(double *arr, long begin, long end, double med_2)
{
  // First sort the orignal array
  double sorted_arr[end-begin];

  long Nt = end - begin;

  for (int i=0; i<Nt; i++) {sorted_arr[i] = arr[begin + i];}

  QuickSort(sorted_arr, 0, Nt - 1); 

  int mid;
  if (Nt % 2 == 0) mid = (int) Nt/2;
  else mid = (int) Nt/2 + 1;
  

  double median = sorted_arr[mid];
  //printf("Sorted median is %f\n",median);
  //printf("Relative error is %f \n", (median-med_2)/median);

  for (int i=0; i<Nt; i++) {arr[begin + i] -= median;}
  return;
}

/*
  Returns the Z (line of sight) and r-ff (orbital) components of the
  lightcurve
  CAUTION: must be called for each sector seperately since folded
  lightcurves may have a different phase. traj_pars must be modified
  accordingly
*/
void Trajectory(double *times, double *traj_pars, double *d_arr, 
  double *Z1_arr, double *Z2_arr, double *rr_arr, 
  double *ff_arr, int Nt)
{
  double M1_cgs = traj_pars[0];
  double M2_cgs = traj_pars[1];
  double P_cgs = traj_pars[2];
  double e = traj_pars[3];
  double inc = traj_pars[4];
  double omega0 = traj_pars[5];
  double T0_cgs = traj_pars[6];

  double Mtot_cgs = M1_cgs + M2_cgs;
  double a = pow(G * Mtot_cgs * SQR(P_cgs) / SQR(2 * PI), 1./3.);

  for (int ind=0; ind<Nt; ind++)
  {
    double t_cgs = times[ind] * SEC_DAY;

    double M = 2.*PI * (t_cgs - T0_cgs) / P_cgs; 
    M = fmod(M,2*PI);
    double EE = M;
    double sin_M = sin(M);

    if(sin_M != 0.0)  EE = M + 0.85*e*sin_M/fabs(sin_M);

    // Kepler's Equation
    for(int j=0; j<7; j++)  EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));

    rr_arr[ind] = a * (1 - e * cos(EE));
    ff_arr[ind] = 2.* atan(sqrt((1. + e)/(1. - e))* tan(EE/2.));

    double cos_omega0_f = cos(omega0 + ff_arr[ind]);
    double sin_omega0_f = sin(omega0 + ff_arr[ind]);
    double cos_inc = cos(inc);
    double sin_inc = sin(inc);

    double ZZ = rr_arr[ind] * sin_omega0_f * sin_inc;
    double d_factor = sqrt(SQR(cos_omega0_f) + SQR(sin_omega0_f * cos_inc));

    d_arr[ind] =  rr_arr[ind] * d_factor;
    Z1_arr[ind] = ZZ * (M2_cgs / Mtot_cgs);
    Z2_arr[ind] = -ZZ * (M1_cgs / Mtot_cgs);
  }
  return;
}

/*
Function to calculate alpha beam
Values taken from Fig 5 (Claret et. al 2020: Doppler beaming factors for white dwarfs
                        main sequence stars and giant stars)
Using values for g=5
Input: log Temperature in Kelvin [NOT NORMALIZED BY SOLAR TEMP]
*/
double Get_Alpha_Beam(double logT)
{
    // Initializing the alpha and temperature values
    double alphas[4] = {6.5, 4.0, 2.5, 1.2};
    double logTs[4] = {3.5, 3.7, 3.9, 4.5};

    // Return endpoints if temperature is outside the domain
    if (logT >= logTs[3]) return 1.2/4;
    if (logT < logTs[0]) return 6.5/4;

    int j = 3;
    while(logT < logTs[j]) j--;

    double alpha_beam = ((alphas[j+1] + (alphas[j+1] - alphas[j]) / (logTs[j+1] - logTs[j]) * (logT - logTs[j+1]))/4 );
    return alpha_beam;
}

/*
Beaming function. Taken from Engel et. al 2020
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'ff' in Trajectory()) (rad)
    alpha_beam  Constant
*/
double Beaming(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double alpha_beam)
{
  double q = M2/M1;
  double fac1 = q / pow(1 + q, 2/3);
  double fac2 = pow(M1, 1./3);
  double fac3 = pow(P, -1./3);
  double fac4 = sin(inc) * cos(omega0 + nu) / sqrt(1 - SQR(e));
  double ppm = 1.e-6;

  double beam = -2830. * alpha_beam * fac1 * fac2 * fac3 * fac4 * ppm;
  return beam;
}


/*
Ellipsoidal variations function. Taken from Engel et. al 2020
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'ff' in Trajectory()) (rad)
    R1:         radius of primary (Rsun)
    a:          semi-major axis (Rsun)
    mu:         limb darkening coefficient
    tau:        gravity darkening coefficient
*/
double Ellipsoidal(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double R1, double a, double mu, double tau)
{
  double alpha_11 = 15 * mu * (2 + tau) / (32 * (3 - mu));
  double alpha_21 = 3 * (15 + mu) * (1 + tau) / (20 * (3 - mu));
  double alpha_2b1 = 15 * (1 - mu) * (3 + tau) / (64 * (3 - mu));
  double alpha_01 = alpha_21 / 9;
  double alpha_0b1 = 3 * alpha_2b1 / 20;
  double alpha_31 = 5 * alpha_11 / 3;
  double alpha_41 = 7*alpha_2b1 / 4;

  double beta = (1 + e * cos(nu)) / (1 - SQR(e));
  double q = M2 / M1;
  // Rotational velocity = Angular velocity at the periapse
  double Prot = P * pow(1 - e, 3./2);

  // "Mean ellipsoidal terms"
  double AM1, AM2, AM3;
  // Sin terms 
  double S1, S3;
  // Cosine terms
  double C2_1, C2_2, C4;

  double ppm = 1.e-6;
  // Flag to include higher order terms in flux
  int enable_const = 1;
  double tot_flux = 0.;

  // First the lowest order Engel Terms
  AM1 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * (1 / SQR(Prot)) * CUBE(R1);
  AM2 = 13435 * 3 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * q / (1 + q) * (1 / SQR(P)) * CUBE(beta*R1);
  C2_1 = 13435 * alpha_21 * SQR(sin(inc)) * (1 / M1) * q / (1 + q) * (1 / SQR(P)) * CUBE(beta * R1) * cos(2*(omega0 + nu));
  
  tot_flux += (AM1 + AM2 + C2_1) * ppm;

  if (enable_const == 1){
      AM3 = 759 * alpha_0b1 * (8 - 40*SQR(sin(inc)) + 35*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3)
      * pow(P, -10./3) * pow(beta * R1, 5);
      S1 = 3194 * alpha_11 * (4 * sin(inc) - 5 * CUBE(sin(inc))) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3)
      * QUAD(beta * R1) * sin(omega0 + nu);
      C2_2 =  759 * alpha_2b1 * (6*SQR(sin(inc)) - 7*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3)
              * (beta * R1) * QUAD(beta * R1) * cos(2*(omega0 + nu));
      S3 = 3194 * alpha_31 * CUBE(sin(inc)) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3) * QUAD(beta * R1) 
              * sin(3*(omega0 + nu));
      C4 = 759 * alpha_41 * QUAD(sin(inc)) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3) * (beta * R1) *
              QUAD(beta * R1) * cos(4*(omega0 + nu));

      tot_flux += (AM3 + S1 + C2_2 + S3 + C4) * ppm;
  }
  return tot_flux; 
}

/*
Reflection function. Taken from Engel et. al (who refer to Faigler & Mazeh 2015)
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'ff' in Trajectory()) (rad)
    R2:         radius of secondary (Rsun)
    alpha_ref1: free paramter for reflection
*/
double Reflection(double P, double M1, double M2, double e, double inc, 
                    double omega0, double nu , double R2, double alpha_ref1)
{
  double q = M2 / M1;
  double beta = (1 + e * cos(nu)) / (1 - SQR(e));
  double ppm = 1.e-6;

  double fac1 = pow(1 + q, -2./3);
  double fac2 = pow(M1, -2./3);
  double fac3 = pow(P, -4./3);
  double fac4 = SQR(beta *R2);
  double fac5 = 0.64 - sin(inc) * sin(omega0 + nu) + 0.18 * SQR(sin(inc)) * (1 - cos(2*(omega0 + nu)));

  return 56514 * alpha_ref1 * fac1 * fac2 * fac3 * fac4 * fac5 * ppm;
}

/*
Median function: Computes the constant terms in the lightcurve model
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    R1:         radius of primary (Rsun)
    R2:         radius of secondary (Rsun)
    a:          semi-major axis (Rsun)
    mu:         limb darkening coefficient
    tau:        gravity darkening coefficient
    alpha_ref1: free paramter for reflection
*/
double Compute_Constant_Terms(double P, double M1, double M2, double e, double inc,
                double omega0, double R1, double R2, double a, 
                double mu, double tau, double alpha_ref1)
{
  double alpha_11 = 15 * mu * (2 + tau) / (32 * (3 - mu));
  double alpha_21 = 3 * (15 + mu) * (1 + tau) / (20 * (3 - mu));
  double alpha_2b1 = 15 * (1 - mu) * (3 + tau) / (64 * (3 - mu));
  double alpha_01 = alpha_21 / 9;
  double alpha_0b1 = 3 * alpha_2b1 / 20;
  double alpha_31 = 5 * alpha_11 / 3;
  double alpha_41 = 7*alpha_2b1 / 4;

  double q = M2 / M1;
  // Rotational velocity = Angular velocity at the periapse
  double Prot = P * pow(1 - e, 3./2);

  // "Mean ellipsoidal terms"
  double AM1, AM2, AM3;


  double ppm = 1.e-6;
  // Flag to include higher order terms in flux
  double median = 0.;

  // First the lowest order Engel Terms
  AM1 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * (1 / SQR(Prot)) * CUBE(R1);
  AM2 = 13435 * 3 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * q / (1 + q) * (1 / SQR(P)) * CUBE(R1 / (1 - SQR(e)));
  AM3 = 759 * alpha_0b1 * (8 - 40*SQR(sin(inc)) + 35*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3)
  * pow(P, -10./3) * pow(R1 / (1 - SQR(e)), 5); 

  double fac1 = pow(1 + q, -2./3);
  double fac2 = pow(M1, -2./3);
  double fac3 = pow(P, -4./3);
  double fac4 = SQR(R2 / (1 - SQR(e)));
  double fac5 = (0.64 + 0.18 * SQR(sin(inc))); 
  double AR = 56514 * fac1 * fac2 * fac3 * fac4 * fac5;

  median += (AM1 + AM2 + AM3 + AR) * ppm;

  return median; 
}

/*
  Extra eclipse contribution with Limb Darkening
  Written by Jeremy Schnittman
*/
double eclipse_area_limb_darkening(double R1, double R2, double d, double u){
  int Nr=50, i;
  double d2, q1, q2, area, dq;
  double q_i, x, phi;
  
  d = fabs(d)/RSUN;
  d2 = d*d;
  q1 = R1*R1;
  q2 = R2*R2;
  dq = q1/Nr;
  if (d > (R1+R2)) area = 0;
  if ((d <= (R1+R2))&(d > 0)) {
    area = 0;
    for (i=0;i<Nr;i++) {
      q_i = (i + 0.5) * dq;
      x = (q_i + d2 - q2) / (2. * sqrt(q_i) * d);
      if (x > 1) phi = 0;
      if (x < -1) phi = PI;
      if ((x >= -1)&(x <= 1)) phi = acos(x);
      area = area + dq * (1. - u * (1. - sqrt(1. - q_i / q1))) * phi;
    }
    area = area/(1.-u/3.);
  }
  if (d == 0) {
    area = 0.;
    for (i=0;i<Nr;i++) {
      q_i = (i+0.5)*dq;
      phi = 0;
      if (q_i <= q2) phi = PI;
      if (q_i > q2) phi = 0;
      area = area + dq*(1.-u*(1.-sqrt(1.-q_i/q1)))*phi;
    }
    area = area/(1.-u/3.);
  }
  return area;
}


/*
Eclipse function - taken from the old likelihood2.c file. Written by Jeremy
Calcualtes the overlapping area of the two stars and fits the radius
Parameters:
    R1:           radius of star 1 (rsun)
    R2:           radius of star 2 (rsun)
    X1:           observed X position of star 1 (cgs)
    X2:           observed X position of star 2 (cgs)
    Y1:           observed Y position of star 1 (cgs)
    Y2:           observed Y position of star 2 (cgs)
IMPORTANT: Make sure to compare Z1 and Z2 before subtracting the relevant flux from the sources
Area returned is in units of solar radius^2. Make sure Xi, Yi, Zi and Ri all have same units!
*/
double Eclipse_Area(double R1, double R2, double d)
{
  double h, r, dc, area, h_sq;
  
  // Function is call by value so values of R1 and R2 aren't swapped in main code
  if (R2 > R1) {
      double temp_ = R1;
      R1 = R2;
      R2 = temp_;
  }
  
  area = 0.;
  d = fabs(d)/RSUN;
  dc = sqrt(R1*R1-R2*R2);
  // Now find the observed overlapping area between the two stars
  if (d >= (R1+R2)) area = 0.;
  if (d < (R1-R2)) area = PI*R2*R2;

  if ((d > dc)&(d < (R1+R2)))
  {
    h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
    h = sqrt(h_sq);

    double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
    double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
    area = Arh1 + Arh2;
  }

  if ((d <= dc)&(d >= (R1-R2)))
  {
    h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
    h = sqrt(h_sq);
    double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
    double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
    area = PI*R2*R2-(-Arh1 + Arh2);
  }
  return area;
}

/*
Function to get the temperature for a star from the log Mass (Msun). I am using tabulated
values given in the TESS portal paper. Note that the final temperature depends on an 
additional scaling parameter and boundary function. Returns temperature in log10K
*/
double _GetT(double logM)
{
  // In solar masses
  double M_nodes[16] = {0.1, 0.26, 0.47, 0.59, 0.69, 0.87,
                        0.98, 1.085, 1.4, 1.65, 2.0, 2.5, 
                        3.0, 4.4, 15., 40.};
  // In log10 K
  double T_nodes[16] = {3.491, 3.531, 3.547, 3.584, 3.644, 3.712,
                        3.745, 3.774, 3.823, 3.863, 3.913, 3.991,
                        4.057, 4.182, 4.477, 4.623};

  double m = pow(10., logM);

  double T;

  // Edge cases
  if (m <= M_nodes[0])
  {
    T = T_nodes[0];
  }
  else if (m >= M_nodes[15])
  {
    T = T_nodes[15];
  }

  // Linear interp otherwise
  else
  {
    for (int j=0; j<16; j++)
    {
      if (m < M_nodes[j])
      {
        T = T_nodes[j-1] + (m - M_nodes[j-1]) * (T_nodes[j] - 
                                T_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1]);
        break;
      }
    }

  }
  return T;
}

/*
Function to get the radius of the star from its mass. Nodes taken from the TESS input catalog
paper and slightly tweaked by John Baker. Note that the final radius depends on an 
additional scaling parameter and boundary function. Returns radius in log10 Rsun
*/
double _GetR(double logM)
{

  // In solar masses
  double M_nodes[10] = {0.07, 0.2, 0.356, 0.655, 0.784, 0.787, 1.377,
                        4.4, 15., 40.};
  // In log10 Rsun
  double logR_nodes[10] = {-0.953, -0.627, -0.423, -0.154, -0.082, -0.087,
                            0.295, 0.477, 0.792, 1.041};

  double m = pow(10., logM);

  double R;

  // Edge cases
  if (m <= M_nodes[0])
  {
    R = logR_nodes[0];
  }

  else if (m >= M_nodes[9])
  {
    R = logR_nodes[9];
  }

  // Linear interp otherwise
  else
  {
    for (int j=0; j<10; j++)
    {
      if (m < M_nodes[j])
      {
        R = logR_nodes[j-1] + (m - M_nodes[j-1]) * (logR_nodes[j] - 
            logR_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1]);
        break;
      }
    }
  }
  return R;
}

/*
Scaling functions for the temperature and the radius. Allows for flexibility in the mass and the
radius. Defined by John Baker; parameters tweaked by Siddhant to make the points lie within
1 sigma
*/
double Envelope_Temp(double logM)
{
  /*The distribution of log10(x/y) ~ N(mu~0, std=0.02264)
  We assume that y(m) = model(m) x 10 ^ (scale x alpha)
  log10(x/y) is a normal distribution therefore we want scale x alpha
  to be a normal distribution. alpha is normally distributed around -1
  and 1 so we just rescale it by multiplying by the std of log10(x/y)
  */
  return 0.0224;
}

double Envelope_Radius(double logM)
{
  double m = pow(10., logM);
  double n = 4.22;
  double slope = 15.68;
  double floor=0.01;
  double corner=1.055;
  double ceil=0.17;

  double boundary = 1/(1/ceil+1/(slope*pow((pow(m, n) + pow(corner, n)), (1/n)) - (slope*corner-floor)));

  return boundary;
}

/*
Compute the radius and Teff for each star
Radii are in solar units
Teffs in K
*/
void Calc_Radii_And_Teffs(double params[],  double *R1, double *R2, double *Teff1, double* Teff2)
{
  double logM1 = params[0];
  double logM2 = params[1];
  double rr1 = params[7];
  double rr2 = params[8];
  double alpha_Teff_1 = params[17];
  double alpha_Teff_2 = params[18];

  *R1 = pow(10., _GetR(logM1) + rr1 * Envelope_Radius(logM1)); 
  *R2 = pow(10., _GetR(logM2) + rr2 * Envelope_Radius(logM2)); 
  *Teff1 = pow(10., _GetT(logM1) + alpha_Teff_1 * Envelope_Temp(logM1));
  *Teff2 = pow(10., _GetT(logM2) + alpha_Teff_2 * Envelope_Temp(logM2));
  return;
}


/*
Full lightcurve calculation
Parameters:
  times:      Time array
  Nt:         Number of points in the time array
  pars:       Parameter array containing the following parameters in order:
      M1:         Mass of primary (log Msun)
      M2:         Mass of secondary (log Msun)
      P:          Period (log days)
      e:          Eccentricity
      inc:        Inclination (rad)
      Omega:      Long of ascending node (rad)
      omega0:     Angle or periastron (rad)
      T0:         Inital Time (days)
      Flux_TESS:  Flux scaling factor
      rr1:        Radius scaling factor of primary (log)
      rr2:        Radius scaling factor of secondary (log)
      mu_1
      tau_1
      mu_2
      tau_2
      alpha_ref_1
      alpha_ref_2
      extra_alpha_beam_1
      extra_alpha_beam_2
      alpha_T1
      alpha_T2
      blending
      flux_tune
  template:   Array to store the lightcurve

*/
void Calculate_Lightcurve(double *times, long Nt, double *pars,
      double *template_)
{
  // Extract the paramters
  double logM1 = pars[0];
  double logM2 = pars[1];
  // Period in seconds
  double P = pow(10., pars[2]) * SEC_DAY;
  double Pdays = pow(10., pars[2]);
  double e = pars[3];
  double cos_inc = pars[4];
  double omega0 = pars[5];
  double T0 = pars[6];  // between 0 and 2 pi
  double rr1 = pars[7];
  double rr2 = pars[8];

  // Beaming coefficients
  int compute_alpha_beam = 1;
  double alpha_beam_1 = 1.;
  double alpha_beam_2 = 1.;

  // Limb and gravity darkening coefficients respectively
  double mu_1 = pars[9];
  double tau_1 = pars[10];
  double mu_2 = pars[11];
  double tau_2 = pars[12];

  // Reflection coefficients
  double alpha_ref_1 = pars[13];
  double alpha_ref_2 = pars[14];

  //extra alphas
  double extra_alpha_beam_1 = exp(pars[15]);
  double extra_alpha_beam_2 = exp(pars[16]);
  double alpha_Teff_1 = pars[17];
  double alpha_Teff_2 = pars[18];

  // lightcurve specific parameters
  double blending = pars[19];
  double flux_tune = pars[20];

  double M1 = pow(10., logM1);
  double M2 = pow(10., logM2);
  double inc = acos(cos_inc);

  // Parameters for the trajectory function
  double M1_cgs = M1 * MSUN;
  double M2_cgs = M2 * MSUN;
  double P_cgs = P;
  double T0_cgs = (T0 / (2 * PI)) * P_cgs;

  double traj_pars[7] = {M1_cgs, M2_cgs, P_cgs, e, inc, omega0, T0_cgs};

  // Compute effective temperature and radius 
  double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

  Calc_Radii_And_Teffs(pars, &R1, &R2, &Teff1, &Teff2);

  // Flux normalization coefficients
  double Norm1, Norm2;
  Norm1 = SQR(R1) * QUAD(Teff1) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));
  Norm2 = SQR(R2) * QUAD(Teff2) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));

  // Set alpha_beam
  if (compute_alpha_beam == 1) 
  {
      alpha_beam_1 = Get_Alpha_Beam(log10(Teff1));
      alpha_beam_2 = Get_Alpha_Beam(log10(Teff2));
  }

  alpha_beam_1 *= extra_alpha_beam_1;
  alpha_beam_2 *= extra_alpha_beam_2;

  // Semi majot axis (cgs calculation)
  double Mtot = (M1+M2)*MSUN;
  double a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  double ar = a / RSUN;

  // Fluxes from the stars stored here
  double Amag1[Nt];
  double Amag2[Nt];

  // Positon arrays for the stars (cylindrical separaion)
  double d_arr[Nt];
  double Z1_arr[Nt];
  double Z2_arr[Nt];

  // Radial separation
  double r_arr[Nt];
  // True anomaly
  double nu_arr[Nt];

  Trajectory(times, traj_pars, d_arr, Z1_arr, Z2_arr, r_arr, nu_arr, Nt);

  // Calculate trajectory and store results in arrays + unitialize fluxes to 0
  for (int i = 0; i<Nt; i++)
  {
    Amag1[i] = 0.;
    Amag2[i] = 0.;

    double beam1 = Beaming(Pdays, M1, M2, e, inc, omega0, nu_arr[i], alpha_beam_1);
    double ellip1 = Ellipsoidal(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R1, ar, mu_1, tau_1);
    double ref1 = Reflection(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R2, alpha_ref_1);

    Amag1[i] = Norm1 * (1 + beam1 + ellip1 + ref1);

    double beam2 = Beaming(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], alpha_beam_2);
    double ellip2 = Ellipsoidal(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R2, ar, mu_2, tau_2);
    double ref2 = Reflection(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R1, alpha_ref_2);

    Amag2[i] = Norm2 * (1 + beam2 + ellip2 + ref2);

    // Eclipse contribution (delta F = F * (ecl area / tot area))
    //double area = Eclipse_Area(R1, R2, d_arr[i]);
    //if (Z2_arr[i] > Z1_arr[i]) Amag2[i] -= area * Norm2 / (PI * SQR(R2));
    //else if (Z2_arr[i] < Z1_arr[i]) Amag1[i] -= area * Norm1 / (PI * SQR(R1));

    // Eclipse with Limb-Darkening
    double area;

	  //if (Z2_arr[i] < Z1_arr[i])
    //{
    //  area = eclipse_area_limb_darkening(R1,R2,d_arr[i],mu_1);
    //  Amag1[i] -= area * Norm1 / (PI * SQR(R1));
	  //}
	  //if (Z2_arr[i] > Z1_arr[i]) 
    //{
	  //  area = eclipse_area_limb_darkening(R2,R1,d_arr[i],mu_2);
	  //  Amag2[i] -= area * Norm2 / (PI * SQR(R2));
	  //} 

    // Full lightcurve
    template_[i] = (Amag1[i] + Amag2[i]);
  }

  // Compare the median
  double median = Norm1 * (1. + Compute_Constant_Terms(Pdays, M1, M2, e, inc, omega0,
                  R1, R2, a, mu_1, tau_1, alpha_ref_1));
  median += Norm2 * (1. + Compute_Constant_Terms(Pdays, M2, M1, e, inc, (omega0+PI),
                  R2, R1, a, mu_2, tau_2, alpha_ref_2));

  // Normalize the lightcurve

  Remove_Median(template_, 0, Nt, median);

  for (int i=0; i<Nt; i++)
  {
    //template[i] -= median;
    template_[i] += 1;
    template_[i] = (1*blending + template_[i]*(1 - blending)) * flux_tune;
  }
  return;
}

/*
Star color calculator
Calculate apparent magnitudes B, G, V, and T given the stellar radii
R1,R2 in cm, T1,T2 in K, and D in pc. 
Sets G, B-V, V-G and G-T mags
*/
void Calc_Mags(double params[],  double D, double *Gmg, double *BminusV, 
              double *VminusG, double *GminusT)
{

  double logM1 = params[0];
  double logM2 = params[1];

  double rr1 = params[7];
  double rr2 = params[8];
  
  double alpha_Teff_1 = params[17];
  double alpha_Teff_2 = params[18];

  double blending = params[19];

  double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

  R1 = pow(10., _GetR(logM1) + rr1 * Envelope_Radius(logM1)); 
  R2 = pow(10., _GetR(logM2) + rr2 * Envelope_Radius(logM2)); 

  Teff1 = pow(10., _GetT(logM1) + alpha_Teff_1 * Envelope_Temp(logM1));
  Teff2 = pow(10., _GetT(logM2) + alpha_Teff_2 * Envelope_Temp(logM2));

  // [Units are Rsun and K respectively]
  // Convert Rsun to cgs
  R1 *= RSUN;
  R2 *= RSUN;

  //wavelength in nm
  double lam[4] = {442,540,673,750};
  double nu[4], f_nu[4];
  double h = 6.626e-27;
  double k = 1.38e-16;
  double pc_cgs = 3.086e18;
  int j;


  for (j=0;j<4;j++)
  {
    nu[j] = C / (lam[j] * 1e-7);
    f_nu[j] = PI * ( R1 * R1 * ( 2. * h * CUBE(nu[j]) / SQR(C) / (exp(h * nu[j] / (k * Teff1)) - 1.)) +
              R2 * R2 * (2. * h * CUBE(nu[j]) / SQR(C) / (exp(h * nu[j] / (k * Teff2)) - 1.)))
              /(SQR(D) * SQR(pc_cgs));

    // Correction from the blending
    f_nu[j] = f_nu[j] / (1 - blending);
  }

  double Bmag = -2.5*log10(f_nu[0])-48.6;
  double Vmag = -2.5*log10(f_nu[1])-48.6;
  double Gmag = -2.5*log10(f_nu[2])-48.6;
  double Tmag = -2.5*log10(f_nu[3])-48.6;

  // Return the differences
  *Gmg = Gmag;
  *BminusV = Bmag - Vmag;
  *VminusG = Vmag - Gmag;
  *GminusT = Gmag - Tmag;
  return;
}

// Function to monitor Roche Lobe overflow (Eggleton 1985)
// Returns Roche Lobe radius over the binary separation
double Eggleton_RL(double q)
{
    return 0.49 * pow(q, 2./3) / (0.6 * pow(q, 2./3) + log(1 + pow(q, 1./3)));
}

// Everything in cgs - I am comparing Roche Lobe to Binary Separation
// at periapse
// Returns 1 if Roche Lobe overflows...
int RocheOverflow(double *pars)
{   
    double M1 = pow(10., pars[0]) * MSUN;
    double M2 = pow(10., pars[1]) * MSUN;
    double q = M1 / M2;
    double period = pow(10., pars[2]) * SEC_DAY;
    double ecc = pars[3];
    double R1 = pow(10., _GetR(pars[0]) + pars[7] * Envelope_Radius(pars[0])) * RSUN; 
    double R2 = pow(10., _GetR(pars[1]) + pars[8] * Envelope_Radius(pars[1])) * RSUN; 
    double Bin_Sep = pow(G * (M1 + M2) * SQR(period) / (4.0*PI*PI), 1./3.);
    // The factor after comes from assuming eccentric orbits
    double RL1_over_sep = Eggleton_RL(q);
    double RL2_over_sep = Eggleton_RL(1/q);


    double R1_over_sep = R1 / (Bin_Sep * (1 - ecc));
    double R2_over_sep = R2 / (Bin_Sep * (1 - ecc));

    int flag = ((RL1_over_sep < R1_over_sep) || (RL2_over_sep < R2_over_sep)) ? 1 : 0;
    return flag;
}

/*
  Finally the likelihood calculation - uses data for multiple
  sectors. Programmed to only use the gmag data
*/
double Log_Likelihood(double all_sector_phases[], double all_sector_fluxes[], 
                      double all_sector_uncertainties[], long int points_per_sector[], 
                      const int NSECTORS, double all_parameters[], 
                      double mag_data[], double mag_err[], const int gmag_flag, 
                      const int color_flag, const int secular_drift_flag)
{
  double logL_net = 0.;

  int skip_samples = 0;
  const int npars_sector = npars_common + npars_unique;

  for (int sector_number = 0; sector_number < NSECTORS; sector_number++)
  {
    double sector_params[npars_sector];
    // Assign common parameters to each sector
    for (int index=0; index<npars_common; index++)
    {
      sector_params[index] = all_parameters[index];
    }

    // Now assign the unique parameters
    for (int index=0; index<npars_unique; index++)
    {
      sector_params[npars_common + index] = all_parameters[npars_common + 
                                            npars_unique * sector_number + 
                                            index];
    }
    // Re-arrange the ordering of the variables if the secular drift flag is on
    if (secular_drift_flag == 1)
    {
      // Current order: logM1, logM2, logP, sigma_r1, sigma_r2, mu_1, tau_1, mu_2, tau_2, alpha_ref_1, alpha_ref_2
      //                alpha_t1, alpha_t2, (e, i, omega, t0, blending, flux_tune, noise_resc)_j
      double __temp[npars_sector];
      __temp[0] = sector_params[0];       __temp[1] = sector_params[1];      __temp[2] = sector_params[2];
      __temp[3] = sector_params[15];      __temp[4] = sector_params[16];     __temp[5] = sector_params[17];
      __temp[6] = sector_params[18];      __temp[7] = sector_params[3];      __temp[8] = sector_params[4];
      __temp[9] = sector_params[5];       __temp[10] = sector_params[6];     __temp[11] = sector_params[7];
      __temp[12] = sector_params[8];      __temp[13] = sector_params[9];     __temp[14] = sector_params[10];
      __temp[15] = sector_params[11];     __temp[16] = sector_params[12];    __temp[17] = sector_params[13];
      __temp[18] = sector_params[14];     __temp[19] = sector_params[19];    __temp[20] = sector_params[20];
      __temp[21] = sector_params[21];

      // And now move back to __temp
      for (int ii=0; ii<npars_sector; ii++)
      {
        //printf("Sector parameters original %d: %f \n", ii, sector_params[ii]);
        sector_params[ii] = __temp[ii];
      }
      for (int ii=0; ii<npars_sector; ii++)
      {
        //printf("Sector parameters changed %d: %f \n", ii, sector_params[ii]);
      }
    }

    // Similary make sector-specific phase, flux and error arrays
    const int Npoints_in_sector = points_per_sector[sector_number];
    double sector_flux[Npoints_in_sector];
    double sector_phase[Npoints_in_sector];
    double sector_uncetainties[Npoints_in_sector];
    double sector_template[Npoints_in_sector];

    for (int index=0; index < Npoints_in_sector; index++)
    {
      sector_flux[index] = all_sector_fluxes[skip_samples + index];
      sector_phase[index] = all_sector_phases[skip_samples + index];
      sector_uncetainties[index] = all_sector_uncertainties[skip_samples + index];
    }

    // Calculate Roche Overflow
    int RocheOverFlowFlag = 0;
    RocheOverFlowFlag = RocheOverflow(sector_params);

    if (RocheOverFlowFlag)
    {
      logL_net = -BIG_NUM;
      //break;
    }

    // Now calculate the lightcurve for the sector
    Calculate_Lightcurve(sector_phase, Npoints_in_sector, sector_params, 
                        sector_template);

    // Finally compute the chi_squared
    double chi2_local = 0.;
    for (int index = 0; index < Npoints_in_sector; index++)
    {
      double point_uncertainty = fmax(1.e-4, sector_uncetainties[index]);
      double residual = (sector_template[index] - sector_flux[index]) /
                        point_uncertainty;
      chi2_local += residual * residual;
    }

    // From the rescaling term
    double ln_noise_resc = sector_params[21];
    chi2_local = chi2_local / exp(2 * ln_noise_resc);
    double logL_resc_term  = - Npoints_in_sector * ln_noise_resc;

    // The the gmag information
    if (gmag_flag)
    {
      double D = mag_data[0];
      double Gmg, BminusV, VminusG, GminusT;

      Calc_Mags(sector_params, D, &Gmg, &BminusV, &VminusG, &GminusT);
      double computed_mags[4] = {Gmg, BminusV, VminusG, GminusT};

      double residual = (Gmg - mag_data[1]) / mag_err[0];
      chi2_local += residual * residual;

      if (color_flag)
      {
        for (int i=0; i<3; i++)
        {
          residual = (computed_mags[i+1] - mag_data[i+2]) / mag_err[i+1];
          chi2_local += residual * residual;
        }
      }
    }
    
      skip_samples += points_per_sector[sector_number];
      logL_net += (-chi2_local/2.0 + logL_resc_term);


  }

  return logL_net;
}

/*
  Return the parameter means and sigmas if they have a gaussian prior
  These must be set in the python script
*/
double Log_Prior(const int NPARS, double *parameter_values, gauss_bounds *gauss_pars)
{
  double LogP = 0.;
  for (int index=0; index<NPARS; index++)
  {
    LogP += log(Gaussian(parameter_values[index], gauss_pars[index].mean, 
    gauss_pars[index].sigma));
  }
  return LogP;
}

/*
  Calculate all the eight components of a lightcurve. These include:
  beaming 1, ellipsoidal 1, reflection 1, reflection 1
  beaming 2, ellipsoidal 2, reflection 2, reflection 2
*/
void Calculate_Lightcurve_Components(double *times, long Nt, double *pars,
      double *beam_arr1, double *ellip_arr1, double *refl_arr1, double *ecl_arr1,
      double *beam_arr2, double *ellip_arr2, double *refl_arr2, double *ecl_arr2)
{
  // Extract the paramters
  double logM1 = pars[0];
  double logM2 = pars[1];
  // Period in seconds
  double P = pow(10., pars[2]) * SEC_DAY;
  double Pdays = pow(10., pars[2]);
  double e = pars[3];
  double cos_inc = pars[4];
  double omega0 = pars[5];
  double T0 = pars[6];  // between 0 and 2 pi
  double rr1 = pars[7];
  double rr2 = pars[8];

  // Beaming coefficients
  int compute_alpha_beam = 1;
  double alpha_beam_1 = 1.;
  double alpha_beam_2 = 1.;

  // Limb and gravity darkening coefficients respectively
  double mu_1 = pars[9];
  double tau_1 = pars[10];
  double mu_2 = pars[11];
  double tau_2 = pars[12];

  // Reflection coefficients
  double alpha_ref_1 = pars[13];
  double alpha_ref_2 = pars[14];

  //extra alphas
  double extra_alpha_beam_1 = exp(pars[15]);
  double extra_alpha_beam_2 = exp(pars[16]);
  double alpha_Teff_1 = pars[17];
  double alpha_Teff_2 = pars[18];

  // lightcurve specific parameters
  double blending = pars[19];
  double flux_tune = pars[20];

  double M1 = pow(10., logM1);
  double M2 = pow(10., logM2);
  double inc = acos(cos_inc);

  // Parameters for the trajectory function
  double M1_cgs = M1 * MSUN;
  double M2_cgs = M2 * MSUN;
  double P_cgs = P;
  double T0_cgs = (T0 / (2 * PI)) * P_cgs;

  double traj_pars[7] = {M1_cgs, M2_cgs, P_cgs, e, inc, omega0, T0_cgs};

  // Compute effective temperature and radius 
  double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

  Calc_Radii_And_Teffs(pars, &R1, &R2, &Teff1, &Teff2);

  // Flux normalization coefficients
  double Norm1, Norm2;
  Norm1 = SQR(R1) * QUAD(Teff1) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));
  Norm2 = SQR(R2) * QUAD(Teff2) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));

  // Set alpha_beam
  if (compute_alpha_beam == 1) 
  {
      alpha_beam_1 = Get_Alpha_Beam(log10(Teff1));
      alpha_beam_2 = Get_Alpha_Beam(log10(Teff2));
  }

  alpha_beam_1 *= extra_alpha_beam_1;
  alpha_beam_2 *= extra_alpha_beam_2;

  // Semi majot axis (cgs calculation)
  double Mtot = (M1+M2)*MSUN;
  double a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  double ar = a / RSUN;

  // Fluxes from the stars stored here
  double Amag1[Nt];
  double Amag2[Nt];

  // Positon arrays for the stars (cylindrical separaion)
  double d_arr[Nt];
  double Z1_arr[Nt];
  double Z2_arr[Nt];

  // Radial separation
  double r_arr[Nt];
  // True anomaly
  double nu_arr[Nt];

  Trajectory(times, traj_pars, d_arr, Z1_arr, Z2_arr, r_arr, nu_arr, Nt);

  // Calculate trajectory and store results in arrays + unitialize fluxes to 0
  for (int i = 0; i<Nt; i++)
  {
    Amag1[i] = 0.;
    Amag2[i] = 0.;

    double beam1 = Beaming(Pdays, M1, M2, e, inc, omega0, nu_arr[i], alpha_beam_1);
    double ellip1 = Ellipsoidal(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R1, ar, mu_1, tau_1);
    double ref1 = Reflection(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R2, alpha_ref_1);

    Amag1[i] = Norm1 * (1 + beam1 + ellip1 + ref1);
    beam_arr1[i] = Norm1 * beam1;
    ellip_arr1[i] = Norm1 * ellip1;
    refl_arr1[i] = Norm1 * ref1;

    double beam2 = Beaming(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], alpha_beam_2);
    double ellip2 = Ellipsoidal(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R2, ar, mu_2, tau_2);
    double ref2 = Reflection(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R1, alpha_ref_2);

    Amag2[i] = Norm2 * (1 + beam2 + ellip2 + ref2);
    beam_arr2[i] = Norm2 * beam2;
    ellip_arr2[i] = Norm2 * ellip2;
    refl_arr2[i] = Norm2 * ref2;

    // Eclipse contribution (delta F = F * (ecl area / tot area))
    double area = Eclipse_Area(R1, R2, d_arr[i]);
    if (Z2_arr[i] > Z1_arr[i])
    {
      Amag2[i] -= area * Norm2 / (PI * SQR(R2));
      ecl_arr2[i] -= area * Norm2 / (PI * SQR(R2));
    }
    else if (Z2_arr[i] < Z1_arr[i])
    {
      Amag1[i] -= area * Norm1 / (PI * SQR(R1));
      ecl_arr1[i] -= area * Norm1 / (PI * SQR(R1));
    }

    // Full lightcurve
    //template[i] = (Amag1[i] + Amag2[i]);
  }

  // Compare the median
  double median = Norm1 * (1. + Compute_Constant_Terms(Pdays, M1, M2, e, inc, omega0,
                  R1, R2, a, mu_1, tau_1, alpha_ref_1));
  median += Norm2 * (1. + Compute_Constant_Terms(Pdays, M2, M1, e, inc, (omega0+PI),
                  R2, R1, a, mu_2, tau_2, alpha_ref_2));

  // Normalize the lightcurve

  Remove_Median(beam_arr1, 0, Nt, median);
  Remove_Median(beam_arr2, 0, Nt, median);
  Remove_Median(ellip_arr1, 0, Nt, median);
  Remove_Median(ellip_arr2, 0, Nt, median);
  Remove_Median(refl_arr1, 0, Nt, median);
  Remove_Median(refl_arr2, 0, Nt, median);
  Remove_Median(ecl_arr1, 0, Nt, median);
  Remove_Median(ecl_arr2, 0, Nt, median);


  for (int i=0; i<Nt; i++)
  {
    beam_arr1[i] += 1;
    beam_arr1[i] = (1*blending + beam_arr1[i]*(1 - blending)) * flux_tune;
    ellip_arr1[i] += 1;
    ellip_arr1[i] = (1*blending + ellip_arr1[i]*(1 - blending)) * flux_tune;
    refl_arr1[i] += 1;
    refl_arr1[i] = (1*blending + refl_arr1[i]*(1 - blending)) * flux_tune;
    ecl_arr1[i] += 1;
    ecl_arr1[i] = (1*blending + ecl_arr1[i]*(1 - blending)) * flux_tune;

    beam_arr2[i] += 1;
    beam_arr2[i] = (1*blending + beam_arr2[i]*(1 - blending)) * flux_tune;
    ellip_arr2[i] += 1;
    ellip_arr2[i] = (1*blending + ellip_arr2[i]*(1 - blending)) * flux_tune;
    refl_arr2[i] += 1;
    refl_arr2[i] = (1*blending + refl_arr2[i]*(1 - blending)) * flux_tune;
    ecl_arr2[i] += 1;
    ecl_arr2[i] = (1*blending + ecl_arr2[i]*(1 - blending)) * flux_tune;
    //template[i] += 1;
    //template[i] = (1*blending + template[i]*(1 - blending)) * flux_tune;
  }
  return;
}