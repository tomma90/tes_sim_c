#include<stdio.h>
#include<math.h>

#define PI 3.14159265
#define RN 1.0
#define TC 0.17
#define ALPHA 100.0
#define RS 0.02
#define L (1.0 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 )
#define IBIAS (33.0 / 10 / 10 / 10 / 10 / 10 / 10)
#define TBASE 0.1                            
#define N 4.0
#define POPT (0.5 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10)  
#define TAU0 (33.0 / 10 / 10 / 10)
#define PSAT (2.5*POPT)
#define G (PSAT * N * TC * TC * TC / (TC * TC * TC *TC - TBASE * TBASE * TBASE * TBASE))
#define C (TAU0 * G)
#define F 0.1   

long double RT(long double T)
{
	long double A = ALPHA * PI * (TC * TC + 1.0) * RN / 2.0 / TC;
	return(RN * ( atan( (T - TC) * A) + PI / 2.0) / PI);
}

long double Pb(long double T, long double Tb)
{
	return(G * (T * T * T * T - Tb * Tb * Tb * Tb) / N / T / T / T);
}

long double V(long double Ib, long double R)
{
	return(Ib * R * RS / (R + RS));
}

long double dIdt(long double t, long double I, long double T, long double R, long double V) 
{
	return((V - I*R - I*RS)/L); 
} 

long double dTdt(long double t, long double I, long double T, long double R, long double Pb, long double Popt) 
{       
	return((I*I*R + Popt - Pb)/C); 
} 

// Finds value of y for a given x using step size h 
// and initial value y0 at x0. 
long double rungeKutta(long double x0, long double yI0, long double yT0, long double x, long double h)
{
	FILE *fileT;
        fileT = fopen("tesDC_T.txt", "w");
	FILE *fileI;
        fileI = fopen("tesDC_I.txt", "w");
	
	// Count number of iterations using step size or 
	// step height h 
	long long int n = (long long int)((x - x0) / h);
	
	long double kI1, kI2, kI3, kI4, kT1, kT2, kT3, kT4; 

	// Iterate for number of iterations 
	long double yI = yI0, yT = yT0;

	long double Popt;
	long double x_0;
	long double Tbase;

	long long int sampling;
	long long int decimate;
	sampling = 1000;
	decimate = (n + 1)/ x / sampling;

	for (long long int i=1; i<=n; i++) 
	{ 
		//Input
		if (x0<=1)
		{
			Popt = POPT;
			x_0 = x0;
		}
		else
		{
			Popt = POPT + 0.01 * POPT * sin(2.0 * PI * F * (x0-x_0));
		}
		Tbase = TBASE; //+ 0.01 * TBASE * sin(2.0 * PI * F * x0);

		// Apply Runge Kutta Formulas to find 
		// next value of y 
		kI1 = h*dIdt( x0, yI, yT, RT(yT), V(IBIAS, RT(yT)) );
		kT1 = h*dTdt( x0, yI, yT, RT(yT), Pb(yT, Tbase), Popt );

		kI2 = h*dIdt( x0 + 0.5*h, yI + 0.5*kI1, yT + 0.5*kT1, RT(yT + 0.5*kT1), V(IBIAS, RT(yT + 0.5*kT1)) );
		kT2 = h*dTdt( x0 + 0.5*h, yI + 0.5*kI1, yT + 0.5*kT1, RT(yT + 0.5*kT1), Pb(yT + 0.5*kT1, Tbase), Popt );

		kI3 = h*dIdt( x0 + 0.5*h, yI + 0.5*kI2, yT + 0.5*kT2, RT(yT + 0.5*kT2), V(IBIAS, RT(yT + 0.5*kT2)) ); 
		kT3 = h*dTdt( x0 + 0.5*h, yI + 0.5*kI2, yT + 0.5*kT2, RT(yT + 0.5*kT2), Pb(yT + 0.5*kT2, Tbase), Popt ); 

		kI4 = h*dIdt( x0 + h, yI + kI3, yT + kT3, RT(yT + kT3), V(IBIAS, RT(yT + kT3)) );
		kT4 = h*dTdt( x0 + h, yI + kI3, yT + kT3, RT(yT + kT3), Pb(yT + kT3, Tbase), Popt ); 

		// Update next value of y 
		yI = yI + (1.0/6.0)*(kI1 + 2*kI2 + 2*kI3 + kI4);
		yT = yT + (1.0/6.0)*(kT1 + 2*kT2 + 2*kT3 + kT4);
		
		if (i%decimate==0)
		{
			long double x0mus = x0*1000000;
			long double yTmK = yT*1000;
			long double yImuA = yI*1000000;
			
			fprintf(fileT, "%Lf\t %Lf\n", x0mus, yTmK);
			fprintf(fileI, "%Lf\t %Lf\n", x0mus, yImuA);

			//printf("%Lf\t%Lf\t%Lf\n", x0mus, yTmK, yImuA);
		}

		// Update next value of x 
		x0 = x0 + h;

	} 

	return x0, yI, yT;

} 

// Driver method 
int main() 
{
	long double t0 = 0, T0 = TC, I0 = V(IBIAS, RT(TC)) / RT(TC), tf = 60, h = 0.000000001;
	long double t, yI, yT = rungeKutta(t0, I0, T0, tf, h);
	printf("\nThe value of I at x is : %Lf\n", yI*1000000);
	printf("\nThe value of T at x is : %Lf\n", yT*1000); 
	return 0; 
} 

