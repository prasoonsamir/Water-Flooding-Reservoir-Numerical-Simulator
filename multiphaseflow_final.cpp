
#include<iostream>
#include<cmath>
#include<string.h>
#include<ctime>
#include<fstream>
#include <vector>

using namespace std;
void print(vector<double>  A);
vector<double> gauss(vector<double> A);
void matrixcoeff(int Nx, int Ny, int Nz, int L, double A[100][100], float Bob, float Cop[50][5][5], float Bw, float Cwp[50][5][5], float Tox[50][5][5], float Twx[50][5][5], float P[50][5][5], string &operating_condition, int welli, int wellj, int wellk, int injwi, int injwj, int injwk, float kro[55][11][10], float krw[55][11][10], float mobo[55][11][10], float mobw[55][11][10], float Viscw, float Viscob, float Ql, float Qw, float WI, float BHP);
void printSw (float Sw[55][10][10], int Ny, int Nx);
void printP (float P[55][10][10], int Ny, int Nx);
void printSwp (float Swproducer[5000], int l);
void transmiss (int Nx, int Ny, int Nz, float P[55][10][10], float mobw[55][10][10], float mobo[55][10][10], float Kx[55][10][10], float Tox[55][10][10], float Twx[55][10][10], float h2, float h3, float h1);

struct date{
int day;
int month;
int year;} datesi, datesf;


 int main()
{  //Read all the input data from input data file by order of top to bottom. Also read the strings in between the data using 'string_read'.....................................................
   float permx_grids, permy_grids, permz_grids, poro_grids, Rs, Pb, den_oil, den_water, den_gas, Ql,Qw,rw,BHP;
    float D[3], permx,permy,permz, phi,topdepth, Prefw, Bw, Cw, Viscw, Pref_rock, Cr, Co,Refdepth, Refpres, WOC, pres, Cp, h1, h2, h3,t,T ;
    int NxNyNz, Nx, Ny, Nz,welli,wellj,wellk,injwi,injwj,injwk,L;
    float RPT[9][4],PVTO[10][3], Swi[9], krwi[9], kroi[9], pcow[9], p[10], Bob,Viscob;
    string string_read,operating_condition;   //Declare the string_read and parameter as string
    ifstream infile;
    infile.open("InputData.txt");//Open the input file named InputData
       infile>>string_read>>string_read>>Nx>>D[0]>>string_read>>Ny>>D[1]>>string_read>>Nz>>D[2];//Reading of dimension data
       NxNyNz=Nx*Ny*Nz;
       infile>>string_read>>WOC;

       infile>>string_read>>topdepth;

      infile>>string_read>>string_read>>permx_grids>>permx>>string_read>>permy_grids>>permy>>string_read>>permz_grids>>permz;

      infile>>string_read>>poro_grids>>phi;

        infile>>string_read>>string_read>>string_read>>string_read>>string_read;
       for (int i=0;i<9;++i)
       {
           for(int j=0;j<4;++j)
           {
               infile>>RPT[i][j];
           }
       }
       for (int i=0;i<9;++i)
       {
           Swi[i]=RPT[i][0];
           krwi[i]=RPT[i][1];
           kroi[i]=RPT[i][2];
           pcow[i]=RPT[i][3];
       }

infile>>string_read>>string_read>>string_read>>string_read>>string_read>>Prefw>>Bw>>Cw>>Viscw;


       infile>>string_read>>string_read>>string_read>>string_read>>Rs>>Pb>>Bob>>Viscob;

       infile>>string_read>>string_read>>string_read>>Pref_rock>>Cr;

       infile>>string_read>>string_read>>string_read>>string_read>>den_oil>>den_water>>den_gas;

       infile>>string_read>>Refdepth>>string_read>>Refpres;

       infile>>string_read>>datesi.month>>datesi.day>>datesi.year;
       infile>>string_read>>datesf.month>>datesf.day>>datesf.year;
       infile>>operating_condition>>Ql;
       infile>>string_read>>BHP;
       infile>>string_read>>Co;
       infile>>string_read>>rw;
       infile>>string_read>>welli>>wellj>>wellk;
       infile>>string_read>>injwi>>injwj>>injwk>>string_read>>Qw;

infile.close();
//Reading of Input data complete...............................................................................................................

//Declare the array and variables................................................................................................................
float P[55][11][10], Pnew[55][11][10], Kx[55][11][10], Ky[55][11][10], Kz[55][11][10], pho[55][11][10], Tox[55][11][10], Toy[55][11][10],Toz[55][11][10],Twx[55][11][10], Twy[55][11][10],Twz[55][11][10],z[10], Swproducer[5000];
float Stable_time,WI,re,V,Pbh,Q[1000], time, Bo[55][11][10],Visco[55][11][10], kro[55][11][10],krw[55][11][10],mobo[55][11][10],mobw[55][11][10],Sn[55][11][10],Cop[55][11][10],Cow[55][11][10],Cwp[55][11][10],Cww,oil_grid, Sw[55][11][10];
//Conversion of input data from field units to SI unit.............................................................................
     pres=Refpres/(1.45*0.0001);
     BHP=BHP/(1.45*0.0001);
     Pb=Pb/(1.45*0.0001);
     Viscob=Viscob/1000;
     Viscw=Viscw/1000;
     den_oil=den_oil*16.0185;
     den_water=den_water*16.0185;
     permx=permx/(1.01*pow(10,15));
     permy=permy/(1.01*pow(10,15));
     permz=permz/(1.01*pow(10,15));
     h1=D[0]/3.28;
     h2=D[1]/3.28;
     h3=D[2]/3.28;
     V=h1*h2*h3;
     rw=rw/3.28;
     Cr=Cr*1.45*pow(10,-4);
     Co=Co*1.45*pow(10,-4);
     Cw=Cw*1.45*pow(10,-4);
     Cp=phi*(Cr+Co);
     Ql=0.00000184*Ql;
     Qw=0.00000184*Qw;
     re=0.14*sqrt(sqrt(permy/permx)*pow(h1,2)+sqrt(permx/permy)*pow(h2,2))/0.5/(pow(permy/permx,0.25)+pow(permx/permy,0.25));
     WI=2*3.14*sqrt(permx*permy)*h3/log(re/rw);                                                                       // Well index
     T=(datesf.year-datesi.year)*365*24*3600+(datesf.month-datesi.month)*30*24*3600+(datesf.day-datesi.day)*24*3600; //Total time for the simulation run
     t=2.4*3600; //Time interval is 1 day
     int n=T/t; // Total number of time steps
     oil_grid=(WOC-topdepth)/D[2];

//Operating condition in the well
if (operating_condition=="PRODUCTIONRATE")
{
  L=Nx*Ny*Nz;
}
else
{//If operating condition is minimum bottomhole pressure
L=Nx*Ny*Nz;
}
//Initialize the pressure and saturation. Also store the permeability values at each direction in a array Kx, Ky and Kz
for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            P[i][j][k]=pres;
            Kx[i][j][k]=permx;
            Ky[i][j][k]=permy;
            Kz[i][j][k]=permz;
            Sw[i][j][k]=Swi[0];
        }
    }
    z[k]=h3*k;
    if (Nz<oil_grid)
    {pres=pres+den_oil*9.8*h3;}
    else
    {
      pres=pres+den_water*9.8*h3;
    }
}
z[Nz]=0;

//Declare the matrix.........
    vector<double> line(L+1,0);
    //vector< vector<double> > A(L,line);
    vector<double> A(L, line);

// Read input data of matrix......
    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++) {
            A[i][j]=0;
        }
    }

//Compute density, formation volume factor, rel perm of water and oil, mobility of oil and water, Cop, Cow, Cwp, Cww for each grid block at every time interval
for (int l=0;l<n; ++l) //Start the time step from l=0 to l=n
{
    for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
           pho[i][j][k]=den_oil*(1+Co*(P[i][j][k]-101325));
           Bo[i][j][k]=Bob/(1+Co*(P[i][j][k]-Pb));
           Sn[i][j][k]=(Sw[i][j][k]-Swi[0])/(1-Swi[0]-0.2);
           krw[i][j][k]=pow(Sn[i][j][k],2);
           kro[i][j][k]=pow(1-Sn[i][j][k],2);
           mobo[i][j][k]=kro[i][j][k]/(Bob*Viscob);
           mobw[i][j][k]=krw[i][j][k]/(Bw*Viscw);
           Cop[i][j][k]=h1*h2*h3/t*(phi*Co/Bob);
           Cow[i][j][k]=-V/t*phi/Bo[i][j][k];
           Cwp[i][j][k]=0;
           Tox[i][j][k]=0;
           Twx[i][j][k]=0;
           Cww=V/t*phi/Bw;

        }
    }
}
int i;
//Calculate transmissibilities of oil and water at the block interfaces in all three dimensions
transmiss(Nx, Ny, Nz, P, mobw, mobo,  Kx, Tox, Twx,  h2,  h3,  h1);


//Setting the matrix coefficient for all grid block
//int *p;
matrixcoeff (Nx, Ny, Nz, L, A, Bob, Cop, Bw, Cwp, Tox, Twx, P, operating_condition, welli, wellj, wellk, injwi, injwj, injwk, kro, krw, mobo, mobw, Viscw, Viscob, Ql, Qw, WI, BHP);
//for (int k=0; k<L; ++k)
//{
//    for (int j=0; j<L; ++j)
//    {
//           int m=j+L*(k);
 //           A[j][k]=p[m];//Pressure solution at new time interval for linearized equation     
 //   }
//}
// Calculate solution of matrix with variable as pressure for different grid blocks using gauss elimination
    vector<double> x(L);
    x = gauss(A);
    
for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
           int m=i+Nx*(j)+Nx*Ny*(k);
            Pnew[i][j][k]=x[m];//Pressure solution at new time interval for linearized equation
        }
    }
}


//Explicit solution of water saturation.............................
 for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {


            Sw[i][j][k]=Sw[i][j][k]+1/Cww*(Twx[i][j][k]*(Pnew[i-1][j][k]-Pnew[i][j][k])+Twx[i+1][j][k]*(Pnew[i+1][j][k]-Pnew[i][j][k]));

    }
}}
//Update of saturation at the production grid block due to production of water
      Sw[welli][wellj][wellk]=Sw[welli][wellj][wellk]+1/Cww*(-Ql*mobw[welli][wellj][wellk]/(mobw[welli][wellj][wellk]+mobo[welli][wellj][wellk]));

     //Update of saturation at the injection grid block due to injection of water
    Sw[injwi][injwj][injwk]=Sw[injwi][injwj][injwk]+1/Cww*Qw;

      for (int k=0; k<Nz; ++k)
{
   for (int j=0; j<Ny; ++j)
   {
       for (int i=0; i<Nx; ++i)
        {
           if (Sw[i][j][k]>0.8)
         {
              Sw[i][j][k]=0.8;
          }
        if (Sw[i][j][k]<0.2)
         {
             Sw[i][j][k]=0.2;
          }

    }
}}

//Subsituting Pold=Pnew and using Pnew for calculating matrix coefficients at new time interval
for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            P[i][j][k]=Pnew[i][j][k];
        }
    }
}
Swproducer[l]=Sw[0][0][0];

}// End of time loop

//Converting pressure in SI unit back to field unit
for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            P[i][j][k]=1.45*0.0001*P[i][j][k];
        }
    }
}

// Call the function print and Print the matrix coefficient in excel file 'linearization','pressure', 'saturation', and 'saturation vs. time' 
    print(A);
    printSw (Sw, Ny, Nx);
    printP (P, Ny, Nx);
    printSwp (Swproducer, l);


return 0;
}







//Calculate transmissibilities of oil and water at the block interfaces in all three dimensions
void transmiss(int Nx, int Ny, int Nz, float P[55][10][10], float mobw[55][10][10], float mobo[55][10][10], float Kx[55][10][10], float Tox[55][10][10], float Twx[55][10][10], float h2, float h3, float h1)
  {
	 for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
           if (i>0)
           {
           if ((P[i][j][k]-P[i-1][j][k])>0)
           {
                mobw[i][j][k]=mobw[i][j][k];
                mobo[i][j][k]=mobo[i][j][k];
        }
            else
            {
               mobw[i][j][k]=mobw[i-1][j][k];
               mobo[i][j][k]=mobo[i-1][j][k];
        }}
          Tox[i][j][k]=2*(h2*h3*Kx[i-1][j][k])*(h2*h3*Kx[i][j][k])/(h1*(h2*h3*Kx[i-1][j][k]+h2*h3*Kx[i][j][k]))*mobo[i][j][k];
          Twx[i][j][k]=2*(h2*h3*Kx[i-1][j][k]*Kx[i][j][k])/(h1*(Kx[i-1][j][k]+Kx[i][j][k]))*mobw[i][j][k];


        }
    }
}}


// Print saturation at end of time loop in excel file 'Saturation'
void printSw (float Sw[55][10][10], int Ny, int Nx)
{
ofstream MyExcelFile("Saturation.csv");
    MyExcelFile << "Result:\t";
for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            MyExcelFile<<Sw[i][j][0]<<",";
        }
        MyExcelFile<<endl;
    }
}
//Print pressure at end of time loop in excel file 'pressure'
void printP (float P[55][10][10], int Ny, int Nx)
{
ofstream MyExcel("pressure.csv");
    MyExcel<< "Result:\t";
for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            MyExcel<<P[i][j][0]<<",";
        }
        MyExcel<<endl;
    }
}
    // Print saturation at producing block vs time in excel file 'Saturation_vs_time'
    void printSwp (float Swproducer[5000], int l)
{  ofstream MyExcel1("saturation_vs_time.csv");
MyExcel1<<"saturation"<<","<<"Time"<<endl;
for (int l=0;l<n;++l)
{
    time=t*(l+1)/3600/24;
    MyExcel1<<Swproducer[l]<<","<<time<<endl;
}}
//Function defined
//void print(vector< vector<double> > A) 
void print(vector<double> A)
{
    int n = A.size();
    ofstream MyExcelFile("linearization.csv");
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            MyExcelFile << A[i][j] << ",";
            if (j == n-1) {
                MyExcelFile << "| ";
            }
        }
        MyExcelFile<<endl;
    }
    cout << endl;
}

// Gauss elimination solver function
//vector<double> gauss(vector< vector<double> > A) 
vector<double> gauss(vector <double> A)
{
    int n = A.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}
//Setting the matrix coefficient for all grid block
void matrixcoeff (int Nx, int Ny, int Nz, int L, vector<double> A[100][100], float Bob, float Cop[50][5][5], float Bw, float Cwp[50][5][5], float Tox[50][5][5], float Twx[50][5][5], float P[50][5][5], string &operating_condition, int welli, int wellj, int wellk, int injwi, int injwj, int injwk, float kro[55][11][10], float krw[55][11][10], float mobo[55][11][10], float mobw[55][11][10], float Viscw, float Viscob, float Ql, float Qw, float WI, float BHP)
{
for (int k=0; k<Nz; ++k)
{
    for (int j=0; j<Ny; ++j)
    {
        for (int i=0; i<Nx; ++i)
        {
            int m=i+Nx*(j)+Nx*Ny*(k);
            A[m][m]=-(Bob*Cop[i][j][k]+Bw*Cwp[i][j][k]+Bob*(Tox[i+1][j][k]+Tox[i][j][k])+Bw*(Twx[i+1][j][k]+Twx[i][j][k]));
            if (m<L-1){A[m+1][m]=Bob*Tox[i+1][j][k]+Bw*Twx[i+1][j][k];}
            if (m<L-1){ A[m][m+1]=Bob*Tox[i+1][j][k]+Bw*Twx[i+1][j][k];}
           // if(m<L-Nx){A[m+Nx][m]=Bo[i][j][k]*Toy[i][j+1][k]+Bw*Twy[i][j+1][k];}
           // if (m<L-Nx){A[m][m+Nx]=Bo[i][j][k]*Toy[i][j+1][k]+Bw*Twy[i][j+1][k];}
           // if (m<L-Nx*Ny){A[m+Nx*Ny][m]=Bo[i][j][k]*Toz[i][j][k+1]+Bw*Twz[i][j][k+1];}
           // if (m<L-Nx*Ny){A[m][m+Nx*Ny]=Bo[i][j][k]*Toz[i][j][k+1]+Bw*Twz[i][j][k+1];}
           A[m][L]=-(Bob*Cop[i][j][k]+Bw*Cwp[i][j][k])*P[i][j][k];
            //+(Bo[i][j][k]*Toz[i][j][k+1]*pho[i][j][k]+Bw*Twz[i][j][k+1]*den_water)*9.8*(z[k+1]-z[k])-(Bo[i][j][k]*Toz[i][j][k]*pho[i][j][k]+Bw*Twz[i][j][k]*den_water)*9.8*(z[k]-z[k-1]);
                 //   Toz[i][j][k+1]*9.8*pho[i][j][k]*(z[k+1]-z[k])-Toz[i][j][k]*9.8*pho[i][j][k]*(z[k]-z[k-1]);


        }
    }
}
if (operating_condition=="PRODUCTIONRATE")
{
int m=welli+Nx*(wellj)+Nx*Ny*(wellk);// grid for production well
int m1=injwi+Nx*(injwj)+Nx*Ny*(injwk);// grid for injection well
A[m][L]=A[m][L]+Ql/(mobw[welli][wellj][wellk]+mobo[welli][wellj][wellk])*(kro[welli][wellj][wellk]/Viscob+krw[welli][wellj][wellk]/Viscw);
cout<<krw[welli][wellj][wellk]<<" "<<kro[welli][wellj][wellk]<<" "<<A[0][L]<<" "<<endl;
A[m1][L]=A[m1][L]-Qw*Bw; // just using the injection rate in the injection well grid with coupling and modifying the matrix coefficient at injection well grid
// Coupling the two equations for production rate and modifying the relevant matrix coefficient
//A[m][m]=A[m][m]-WI*(mobo[welli][wellj][wellk]*Bo[welli][wellj][wellk]+mobw[welli][wellj][wellk]*Bw);
//A[m][L-1]=WI*(mobo[welli][wellj][wellk]*Bo[welli][wellj][wellk]+mobw[welli][wellj][wellk]*Bw);
//A[L-1][m]=-WI*(mobo[welli][wellj][wellk]+mobw[welli][wellj][wellk]);
//A[L-1][L-1]=WI*(mobo[welli][wellj][wellk]+mobw[welli][wellj][wellk]);
//A[L-1][L]=-Ql;

//A[m][m]=A[m][m]-WI*(mobw[injwi][injwj][injwk]*Bw);
//A[m][L-1]=WI*(mobw[injwi][injwj][injwk]*Bw);
//A[L-1][m]=-WI*(mobw[injwi][injwj][injwk]);
//A[L-1][L-1]=WI*(mobo[injwi][injwj][injwk]+mobw[injwi][injwj][injwk]);
//A[L-1][L]=Qw;
}
else
{//If operating condition is minimum bottomhole pressure and modifying the relevant matrix coefficient
  int m=welli+Nx*(wellj)+Nx*Ny*(wellk);
A[m][m]=A[m][m]-WI*(mobo[welli][wellj][wellk]*Bob+mobw[welli][wellj][wellk]*Bw);
A[m][L]=A[m][L]- WI*(mobo[welli][wellj][wellk]*Bob+mobw[welli][wellj][wellk]*Bw)*BHP;
}
}
