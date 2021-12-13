#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#define PI 3.14159265
#define RN 1.0
#define TC 0.17
#define ALPHA 100.0
#define RS 0.02
#define L (1.0 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 )
#define IBIAS (33.0 / 10 / 10 / 10 / 10 / 10 / 10) * pow(2, 0.5)
#define TBASE 0.1                            
#define N 4.0
#define POPT (0.5 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10 / 10)  
#define TAU0 (33.0 / 10 / 10 / 10)
#define PSAT (2.5*POPT)
#define G (PSAT * N * TC * TC * TC / (TC * TC * TC *TC - TBASE * TBASE * TBASE * TBASE))
#define C (TAU0 * G)
#define F 2.0
#define FTONE 1000000.0
#define FLC 1000000.0
#define LLC (65.0 / 10 / 10 / 10 / 10 / 10 / 10)
#define CLC (1 / LLC / (2.0 * PI * FLC) / (2.0 * PI * FLC))   

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

long double dJdt(long double t, long double J, long double I, long double T, long double R, long double V)
{
	long double temp = I/CLC;
	return((V-J*(R)-temp)/(L+LLC));
}

long double dIdt(long double t, long double J, long double I) 
{
	return(J); 
} 

long double dTdt(long double t, long double I, long double T, long double R, long double Pb, long double Popt) 
{       
	return((I*I*R + Popt - Pb)/C); 
} 

// Finds value of y for a given x using step size h 
// and initial value y0 at x0. 
long double rungeKutta(long double x0, long double yJ0, long double yI0, long double yT0, long double x, long double h)
{	
	FILE *fileT;
	FILE *fileI;
	
	// Count number of iterations using step size or 
	// step height h 
	long long int n = (long long int)((x - x0) / h);
	
	long double kI1, kI2, kI3, kI4, kT1, kT2, kT3, kT4, kJ1, kJ2, kJ3, kJ4; 

	// Iterate for number of iterations 
	long double yJ = yJ0, yI = yI0, yT = yT0;

	long long int i_new = 1;
	int j = 1;

	for (long long int i=1; i<=n; i++) 
	{ 
		//Input
		long double Ibase = IBIAS * 2.0 * PI * FLC * cos(2.0 * PI * FLC * x0);
		long double Popt = POPT; //+ 0.01 * POPT * sin(2.0 * PI * F * x0);
		long double Tbase = TBASE; //+ 0.01 * TBASE * sin(2.0 * PI * F * x0);

		// Apply Runge Kutta Formulas to find 
		// next value of y
		kI1 = h*dIdt( x0, yJ, yI ); 
		kJ1 = h*dJdt( x0, yJ, yI, yT, RT(yT), V(Ibase, RT(yT)) );
		kT1 = h*dTdt( x0, yI, yT, RT(yT), Pb(yT, Tbase), Popt );

		kI2 = h*dIdt( x0 + 0.5*h, yJ + 0.5*kJ1, yI + 0.5*kI1);
		kJ2 = h*dJdt( x0 + 0.5*h, yJ + 0.5*kJ1, yI + 0.5*kI1, yT + 0.5*kT1, RT(yT + 0.5*kT1), V(Ibase, RT(yT + 0.5*kT1)) );
		kT2 = h*dTdt( x0 + 0.5*h, yI + 0.5*kI1, yT + 0.5*kT1, RT(yT + 0.5*kT1), Pb(yT + 0.5*kT1, Tbase), Popt );
		
		kI3 = h*dIdt( x0 + 0.5*h, yJ + 0.5*kJ2, yI + 0.5*kI2);
		kJ3 = h*dJdt( x0 + 0.5*h, yJ + 0.5*kJ2, yI + 0.5*kI2, yT + 0.5*kT2, RT(yT + 0.5*kT2), V(Ibase, RT(yT + 0.5*kT2)) ); 
		kT3 = h*dTdt( x0 + 0.5*h, yI + 0.5*kI2, yT + 0.5*kT2, RT(yT + 0.5*kT2), Pb(yT + 0.5*kT2, Tbase), Popt ); 

		kI4 = h*dIdt( x0 + h, yJ + kJ3, yI + kI3);
		kJ4 = h*dJdt( x0 + h, yJ + kJ3, yI + kI3, yT + kT3, RT(yT + kT3), V(Ibase, RT(yT + kT3)) );
		kT4 = h*dTdt( x0 + h, yI + kI3, yT + kT3, RT(yT + kT3), Pb(yT + kT3, Tbase), Popt ); 

		// Update next value of y 
		yI = yI + (1.0/6.0)*(kI1 + 2*kI2 + 2*kI3 + kI4);
		yJ = yJ + (1.0/6.0)*(kJ1 + 2*kJ2 + 2*kJ3 + kJ4);
		yT = yT + (1.0/6.0)*(kT1 + 2*kT2 + 2*kT3 + kT4);
		
		long long int sampling = 25000000;
		long long int decimate = (n + 1) / sampling;
                
		char fT[1000];
		char num[1000];
		sprintf(num, "%d", j);
		char txt[1000];
		strcpy(txt, ".txt");
		strcpy(fT, "tesAC_T_");
		strcat(fT, num);
		strcat(fT, txt);
		char fI[1000];
		strcpy(fI, "tesAC_I_");
		strcat(fI, num);
		strcat(fI, txt);

		if (i%decimate==0)
		{
			long double yTmK = yT*1000;
			long double yImuA = yI*1000000;


			if (i_new > (j-1)*10000 && i_new <= j*10000)
			{	
       				fileT = fopen(fT, "a");
				fileI = fopen(fI, "a");
			
				fprintf(fileT, "%Lf\t %Lf\n", x0, yTmK);
				fprintf(fileI, "%Lf\t %Lf\n", x0, yImuA);

				fclose(fileI);
				fclose(fileT);
			}
			

			if (i_new==j*10000)
			{
				j = j + 1;
			}
	
			i_new = i_new + 1;
		}		

		// Update next value of x 
		x0 = x0 + h;

	} 

	return x0, yI, yT;

} 

// Driver method 
int main() 
{
	long double t0 = 0, T0 = TC, J0 = 0, I0 = V(IBIAS, RT(TC)) / RT(TC), tf = 1, h = 0.000000001;
	long double t, yI, yT = rungeKutta(t0, J0, I0, T0, tf, h);
	printf("\nThe value of I at x is : %Lf\n", yI*1000000);
	printf("\nThe value of T at x is : %Lf\n", yT*1000); 
	return 0; 
} 

