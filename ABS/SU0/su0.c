

/* This program generates the model that goes to CPLEX with consideration of spectrum efficiency*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include "header/SpectrumEfficiency.h"
//#include "header/PowerSaving.h"
#include "header/BPSK.h"
#include "header/64QAM.h"
#include "header/16QAM.h"
#include "header/RANDOM.h"
#include "header/Macro.h"
#include "header/poisson.h"
#include "header/ONUserSeEOptim.h"
#include "header/ComparisonModel.h"
#include "header/comparison1.h"
#include "header/comparison2.h"
#include "header/comparison4.h"
#include "header/comparison5.h"
#include "header/OptModel.h"
#include "header/NoMod.h"
#include "header/EEMax.h"
#include "header/SEMax.h"
#include "header/SpecParti.h"
#include "header/SmallUser.h"
#include "header/association.h"
#include "header/zooming.h"
#include "header/switch.h"
#include "header/MCS.h"

#include "header/OptModelABS.h"
#include "header/OptFixedABSpattern.h"
#include "header/FixedABSpattern.h"
#include "header/locOpt.h"
#include "header/SleepingStrategies.h"

#define map_x 3000 //x-axis of map
#define map_y 3000	//y-axis of map
#define Ntp 360	//number of TPs
#define lambda 30	//mean arrival users per cell
#define space 20	//minimum unit on map, should be 20m
#define Rbs 300	//radius of BS, should be 5000km
#define rbs 30	//the radius of a micro BS in meter
#define MNbs 38	//assumed maximum number of BSs including micro and macro BSs
#define MBR 990 * 2	//maximum bandwidth requirement QPSK 3/4
#define mNbs 10	//number of micro BSs
//table:parameters per channel bandwidth
#define BW 1428571	//10 MHz channel bandwidth
//#define DSt 1680000	//total data subcarriers for DL subframes within 1ms (7 DL symbols per subframe) (5040*1000 in our TMC paper)
#define DSt 7142	//total resource blocks in 1000 DL subframe within 1 second 
#define PS 120	//pilot subcarriers
#define FFT 341	//FFT size
#define N 1.12	//n in thermal noise equation
//inputs
#define NF 6	//Rx Noise Figure
#define BH 30	//BS antenna height
#define mBH 5	//micro BS antenna height
#define MH 1	//MS antenna height
#define Gta 16	//BS antenna gain
#define mGta 2	//micro BS antenna gain
#define Gto 9	//other DL Tx gain
#define mGto 0	//other DL Tx gain of micro BSs
#define Gra 2	//Rx antenna gain
#define Gro 0	//Rx other gain
#define FFM 3	//fast fading margin

#define MP 6.3		//UNIT 1mW=30dBm maximum power of a BS: 20W
#define SBP 75.00	//Power for BS in sleep mode
#define ABP 130.0	//Basic power for BS in active mode
#define BP 0.0	//UNIT(%) BLOCKING PROBABILITY
#define mP 0.1	//maximum transmission power of micro BS: 10W
#define mSBP 2.8	//Power for micro BSs in sleep mode
#define mABP 4.8	//Basic power for micro BSs in active mode

//area of urban and suburban
#define X1urban 0
#define X2urban 150
#define	Y1urban 0
#define Y2urban 150
#define X1suburban 0
#define X2suburban 150
#define Y1suburban 0
#define Y2suburban 150
#define UrbanUser 0.075
#define SubUser 0.125
#define RuralUser 0.8

#define LOSmin 20	//los distance
#define NLOSmin 1000	//nlos disance

#define Modulation 3	//number of modulation schemes
#define MaxCellSize 420.0	//the maximum macro cell size
#define maxCellSize 120.0	//the maximum micro cell size

#define limit 0.429

double UEacrossArea = 7.2;
double UEDensity = 50.0;
double gridSize = 10.0;
int UEhotSpotGrids = 4;
//double hotSpotUEDensity = UEDensity * 2;
int Xbs[MNbs]; int Ybs[MNbs]; // x- and y-coordinates of BS sites
int Xtp[Ntp]; int Ytp[Ntp]; // x- and y-coordinates or TPs
double D[MNbs][Ntp];	//distance from BS to TP
double* ptrD;	//pointer to D
double BR[Ntp];	//bandwidth requirement of TPs
double* ptrBR;	//pointer to BR
int Mm[MNbs]={0};	//micro BSs number in indivisual macro BSs
int* ptrMm;
//table:parameters per modulation scheme
int AdjacentBSs[MNbs][Ntp];	//adjacent BSs for users
double Sensitivity[Ntp][Modulation];	//Sensitivity for TP
double DP[Ntp][Modulation];	//Data Bit per Symbol for TP
int minM[MNbs-mNbs][MNbs];	//micro BS in Macro BS
int* ptrAdjacentBSs;	//pointer to AdjacentBSs[][]
double* ptrSensitivity;	//pointer to Sensitivity[]
double* ptrDP;	//pointer to DP[]
int* ptrminM;	//pointer to minM[][]

//table:urban corrections
int LT[Ntp];	//location type of TP
int* ptrLT;	//pointer to LT[]
int BPL[Ntp];	//building penetration loss to TP
int* ptrBPL;	//pointer to BPL

int DS[Ntp][Modulation];	//number of data subcarriers by TP
int* ptrDS;	//pointer to DS[][]

double PL[MNbs][Ntp];	//pathloss between BS and TP
int PathLossType[MNbs][Ntp];	//pathloss type of TP
double* ptrPL;	//pointer to PL[][]
int* ptrPLT;	//pointer to PLT[]

double P[MNbs][Ntp][Modulation];	//power for the transmission from BS to MS
double* ptrP;
 //double Mbt[_t];

void checkAdjacentBSs(int* Xbs, int* Ybs, int* Xtp, int* Ytp, int Tbs);
void specifyTypes(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void bounds(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void checkMSsites(int Tbs);
void FixedCellSize(int Tbs, int* DS, double* ptrD, double* P, double* ptrBR, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int* ptrminM, int nUser);
void power(int Tbs, double TN);
void fillCoordinatesBSs(int* Xbs, int* Ybs, int Gbs_x, int Gbs_y, int Nbs_x1, int Nbs_x2, int Nbs_y1, int Nbs_y2, int length, int Tbs, int* ptrMm, int* ptrminM, int Px, int Py);
void fillCoordinatesTPs(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs, int nUser);
double distance(double xa, double ya, double xb, double yb);
void fillDistance(int* Xbs, int* Ybs, int* Xtp, int* Ytp, double* ptrD, int Tbs);
//void fillMs(double* ptrBR);
double ThermalNoise();
void modulation(double TN, double* ptrSensitivity, double* ptrDP);//double* ptrBR, ,int* ptrModulation,  double* ptrSNR);
void LocationType(int* ptrLT, int* Xtp, int* Ytp);
void bpl(int* ptrLT, int* ptrBPL);
void DataSubcarriers(int* ptrDS, double* ptrBR, double* ptrDP);
void PathLoss(double* ptrD, int* ptrPLT, int* ptrLT, int Tbs, double* ptrPL);
void UECoordinates(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs, int nUser);


void constraint9(int Tbs, FILE* NoExistingUser, FILE* Feasibility);
void constraint8(int Tbs, FILE* NoExistingUser);
void constraint7(int Tbs, FILE* NoExistingUser);
void constraint6(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint5(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint4(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint3(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, int nUser, FILE* Feasibility);
void constraint2(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint1(int Tbs, int* ptrDS,FILE* PowerSavinCPLEX,int* ptrMm, int* ptrminM, FILE* NoExistingUser, FILE* Feasibility);
void objective(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint10(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility);
void constraint11(int Tbs,FILE* PowerSavingCPLEX,int nUser, FILE* Feasibility);
//-------------------------------------------------------
main() {
 int Gbs_x=map_x/Rbs; // Number of grids on y=1
 int Gbs_y=map_y/Rbs; // Number of grids on x=1
 int Nbs_x1=floor((Gbs_x-2)/3)+1;	//Number of BSs on y=1
 int Nbs_x2=floor(((map_x-(Rbs/2))/Rbs)/3);	//Number of BSs on y=2
 int Nbs_y1=floor(Gbs_y/2);	//Number of BSs on x=1
 int Nbs_y2=floor(((map_y-Rbs)/Rbs)/2);	//Number of BSs on x=2
	if(Nbs_x1<1||Nbs_x2<1||Nbs_y1<1||Nbs_y2<1){
		printf("Initial Map Error. Program Stops!");
		exit(1);
	}
 int Tbs=(Nbs_x1*Nbs_y1)+(Nbs_x2*Nbs_y2)+mNbs;	//Total BS on map
 int length=Rbs/space;	//real length
 int Px=map_x/space;	//number of points on x-axis
 int Py=map_y/space;	//number of points on y-axis
 double TN;	//Thermal Noise
 int i,j,k;
//printf("%lf",(float)736/800);
//exit(1);
 ptrD=&D[0][0];
 ptrBR=&BR[0];
	ptrMm=&Mm[0];
 ptrSensitivity=&Sensitivity[0][0];
 ptrDP=&DP[0][0];
 ptrAdjacentBSs=&AdjacentBSs[0][0];
 ptrBPL=&BPL[0];
 ptrLT=&LT[0];
 ptrDS=&DS[0][0];
 ptrPLT=&PathLossType[0][0];
 ptrPL=&PL[0][0];
 ptrP=&P[0][0][0];
	ptrminM=&minM[0][0];
//printf("\n%d\n",Nbs_x1);printf("\n%d\n",Nbs_x2);printf("\n%d\n",Nbs_y1);printf("\n%d\n",Nbs_y2);printf("\n%d\n",Tbs);
  
 //for(counter=0;counter<_t;counter++)
	 //totalr+=R[counter];//totalr = 4+6+8+5+2
	int maxAUser=Ntp/(Tbs-mNbs);	//maximum arrival users per cell
	int nUser;	//total number of new coming users
 //checkInputs();
 eNodeBCoordinates(Xbs, Ybs, &Tbs);
 //fillCoordinatesBSs(Xbs,Ybs,Gbs_x,Gbs_y,Nbs_x1,Nbs_x2,Nbs_y1,Nbs_y2,length,Tbs,ptrMm,ptrminM, Px, Py); 
 
 /*for(i=1;i<=Tbs;i++)
	 printf("Xbs[%d],Ybs[%d]=(%d,%d)\n",i,i,Xbs[i],Ybs[i]);
*/nUser = Ntp;//(Tbs-mNbs)*PoissonRandomNumber(lambda,maxAUser);
 UECoordinates(Xtp, Ytp, Px, Py, Tbs, Xbs, Ybs, nUser);
 //fillCoordinatesTPs(Xtp,Ytp,Px,Py,Tbs,Xbs,Ybs,nUser);
//  SmallUser(Xtp,Ytp,Px,Py,Tbs,Xbs,Ybs,nUser,mNbs,space,maxCellSize,Rbs,Ntp);

 checkAdjacentBSs(Xbs,Ybs,Xtp,Ytp,Tbs);
 fillDistance(Xbs,Ybs,Xtp,Ytp,ptrD,Tbs);/*
 printf("\n");
for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++)
		printf("D[%d][%d]=%.2lf\n",i,j,D[i][j]);
printf("\n");	*/
 fillMs(ptrBR, Xtp, Ytp);
/* for(i=1;i<=Ntp;i++)
	 printf("BR[%d]:%.2lf Mbps\n",i,BR[i]);
 printf("\n");*/
 TN=ThermalNoise();
 //printf("%.2lf\n",TN);
 //printf("\n");
 modulation(TN,ptrSensitivity,ptrDP);//,ptrBR,Modulation,SNR);
 //for(i=1;i<=Ntp;i++)
//	 printf("TP[%d] with modulation %d, DP: %.2lf, SNR: %.2lf, S: %.2lf\n",i,Modulation[i],DP[i],SNR[i],Sensitivity[i]);
 //printf("\n");
 //LocationType(ptrLT,Xtp,Ytp);
 bpl(ptrLT,ptrBPL);
 //for(i=1;i<=Ntp;i++)
//	 printf("TP[%d]'s location type is %d, bpl is %d\n",i,LT[i],BPL[i]);
 DataSubcarriers(ptrDS,ptrBR,ptrDP);
 //printf("\n");

 	
	//for(i=1;i<=Ntp;i++)
	// printf("The number of new coming users is %d\n",nUser);
	 //system("pause");
 PathLoss(ptrD,ptrPLT,ptrLT,Tbs,ptrPL);/*
  for(i=0;i<Tbs;i++)
    for(j=0;j<Ntp;j++)
	printf("PL[%d][%d}=%g\n",i,j,PL[i][j]);*/
power(Tbs,TN);
FixedCellSize(Tbs,ptrDS,ptrD,ptrP,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrminM,nUser);
/*
FILE *PowerSavingCPLEX;
if((PowerSavingCPLEX=fopen("PowerSavingCPLEX","w"))==NULL)
	printf("\nerror!Fail to open file!");
else
	printf("\nOpen PowerSavingCPLEX successfully!\n");
fprintf(PowerSavingCPLEX,"This is the input to CPLEX for power saving model.\n");
	FILE *NoExistingUser;
	if((NoExistingUser=fopen("NoExistingUser","w"))==NULL)
		printf("\nerror!Fail to open the file of NoExistingUser.txt");
	else
		printf("\nOpen NoExistingUser.txt successfully!\n");
	fprintf(NoExistingUser,"This is the input to CPLEX for no exisinting users.\n");
	
	FILE *Feasibility;
	if((Feasibility=fopen("Feasibility","w"))==NULL)
		printf("\nerror!Fail to open the file of Feasibility.txt");
	else
		printf("\nOpen Feasibility.txt successfully!\n");
	fprintf(Feasibility,"This is the input to CPLEX for feasibility check.\n");
	
 	objective(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	
 	fprintf(NoExistingUser,"\nst\n");
fprintf(PowerSavingCPLEX,"\nst\n");
 printf("\nst\n");
 
	constraint1(Tbs,ptrDS,PowerSavingCPLEX,ptrMm,ptrminM,NoExistingUser,Feasibility);
	constraint2(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	constraint3(Tbs,PowerSavingCPLEX,NoExistingUser,nUser,Feasibility);
	constraint4(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	constraint5(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	constraint10(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	constraint6(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
 constraint7(Tbs,NoExistingUser);
 constraint8(Tbs,NoExistingUser);
	constraint9(Tbs,NoExistingUser,Feasibility);
	constraint11(Tbs,PowerSavingCPLEX,nUser,Feasibility);
	fprintf(NoExistingUser,"bounds\n");
fprintf(PowerSavingCPLEX,"bounds\n");
 printf("bounds\n");
	bounds(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
	double ActivePower=ABP;
	double SleepPower=SBP;
	specifyTypes(Tbs,PowerSavingCPLEX,NoExistingUser,Feasibility);
fprintf(PowerSavingCPLEX,"end\n");
 printf("end\n");
 	fprintf(NoExistingUser,"end\n");
fclose(PowerSavingCPLEX);
	fclose(NoExistingUser);
	fclose(Feasibility);*/
 checkMSsites(Tbs);
	//heuristic(Tbs,ptrP,Ntp,MP,DSt,BP,ptrDS,Xbs,Ybs,Xtp,Ytp,ptrD,BR,Modulation,mABP,mSBP);
	//heuristicXUser(Tbs,ptrP,Ntp,MP,DSt,BP,ptrDS,Xbs,Ybs,Xtp,Ytp,ptrD,BR,Modulation,mABP,mSBP,mNbs,mP,ptrAdjacentBSs,ptrminM,ptrMm,nUser,ABP,SBP);
	printf("===================================================================================================================================\n");
	Sheuristic(Tbs,ptrP,Ntp,MP,DSt,BP,ptrDS,Xbs,Ybs,Xtp,Ytp,ptrD,BR,Modulation,ABP,SBP,mABP,mSBP,mNbs,mP,ptrAdjacentBSs,ptrminM,ptrMm,nUser); 

FILE *outfile, *outfile1;
 if ((outfile=fopen("outfile1.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen file successfully!\n");
fprintf(outfile,"BS%dMS%dBP%g\n",Tbs,Ntp,BP);
for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
for(k=0;k<Modulation;k++)
	fprintf(outfile,"%d. BS[%d]=%d,%d MS[%d]=%d,%d with modulation %d DS=%d power=%g mW BW=%g Mbps \n",i*Ntp*Modulation+j*Modulation+k,i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],k,DS[j][k],P[i][j][k],BR[j]);
fclose(outfile);

if ((outfile1=fopen("coordinates.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen coordinates.txt successfully!\n");

fprintf(outfile1,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
//for(i=0;i<Tbs;i++){
for(j=0;j<Ntp;j++)
//for(k=0;k<Modulation;k++)
//fprintf(outfile1,"%d. BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BW %g Mbps with modulation %d\n",i*Ntp*Modulation+j*Modulation+k,i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],DS[j][k],P[i][j][k],BR[j],k);
fprintf(outfile1,"MS[%d] %d %d\n",j,Xtp[j],Ytp[j]);
//}
fprintf(outfile1,"\n\n");
  for (i=0;i<Tbs;i++)
    if (i < Tbs - mNbs){
      fprintf(outfile1,"\"Macro-%d\" %d %d %d %d\n", i, Xbs[i], Ybs[i], Rbs / space, (int)(MaxCellSize / space));
      if (i == Tbs - mNbs - 1)
	fprintf(outfile1, "\n\n");
    }
    else
      fprintf(outfile1,"\"Femto-%d\" %d %d %d %d\n", i - Tbs + mNbs, Xbs[i], Ybs[i], rbs / space, (int)(maxCellSize / space));
  fprintf(outfile1,"\n\nplot \"coordinates.txt\" index 1:1 using 2:3 title \"Macros\" with points, \"coordinates.txt\" index 2:2 using 2:3 title \"Femtos\" with points, \"coordinates.txt\" index 1:1 using 2:3:4 title \"Macros\" with circles, \"coordinates.txt\" index 2:2 using 2:3:4 title \"Femtos\" with circles, \"coordinates.txt\" index 0:0 using 2:3 title \"UEs\" with points");
fclose(outfile1);
	FILE *parameter;
	if((parameter=fopen("parameter.txt","w"))==NULL){
		printf("\n\nerror! Fail to open file!");
		exit(1);
	}
	else
		printf("\n\nOpen parameter.txt successfully!\n"); //printf("The number of new coming users is %d\n",nUser);
	fprintf(parameter,"Tbs=%d\nNtp=%d\nBP=%g\nDSt=%d\nmNbs=%d\nModulation=%d\nmP=%lf\nMP=%lf\nSBP=%lf\nmSBP=%lf\nmABP=%lf\nABP=%lf\n",Tbs,Ntp,BP,DSt,mNbs,Modulation,mP,MP,SBP,mSBP,mABP,ABP);
	fprintf(parameter, "The number of new coming users is %d\n",nUser);
	for(i=0;i<Tbs-mNbs;i++)
		for(j=0;j<Tbs;j++)
		  if (minM[i][j] == 1)
			fprintf(parameter,"minM[%d][%d]=%d\n",i,j,minM[i][j]);
	for(i=0;i<Tbs;i++){
	  int UEInCoverage = 0;
		for(j=0;j<Ntp;j++)
		  if (AdjacentBSs[i][j] == 1){
			fprintf(parameter,"AdjacentBSs[%d][%d]=%d\n",i,j,AdjacentBSs[i][j]);		
			UEInCoverage ++;
		  }
	  if (i >= Tbs- mNbs && UEInCoverage > 0){
	    fprintf(parameter, "UEs in coverage of eNodeB%d: %d\n", i, UEInCoverage);
	    callingTimesCalculation(UEInCoverage, parameter);
	  }
	}
	fclose(parameter);
	
	localOptimal(Tbs, nUser, Ntp, Modulation, mNbs, DSt, mSBP, mABP, ABP, SBP, MP, mP, BP);
 return;
}

eNodeBCoordinates(int* Xbs, int* Ybs, int* Tbs){
  FILE* BSCoIn;
  if ((BSCoIn = fopen("eNodeB.txt", "r")) == NULL){
    printf("Fail to open file\n");
    exit(1);
  }
  else 
    printf("Open eNodeB.txt successfully!\n");
  char id[1000];
  int label, X, Y;
  int macroNum = 0;
  while (fscanf(BSCoIn, "%s", &id) != EOF){
    
    fscanf(BSCoIn, " \"Macro-%d\" %d %d", &label, &X, &Y);
    *(Xbs + label) = X * 3.7 / (double)space;
    *(Ybs + label) = (632 - Y) * 3.7 / (double)space;
   // printf("%d %d %d\n", label, X, Y);
  }
  fclose(BSCoIn);
  macroNum = label + 1;
  FILE* FeCoIn;
  if ((FeCoIn = fopen("femto.txt", "r")) == NULL){
    printf("Fail to open file\n");
    exit(1);
  }
  else 
    printf("Open femto.txt successfully!\n");

  while (fscanf(FeCoIn, "%s", &id) != EOF){
    
    fscanf(FeCoIn, " \"Femto-%d\" %d %d", &label, &X, &Y);
    *(Xbs + label + macroNum) = X * 3.7 / (double)space;
    *(Ybs + label + macroNum) = (632 - Y) * 3.7 / (double)space;
    //printf("%d %d %d\n", label, X, Y);
  }
  fclose(FeCoIn);
  *Tbs = macroNum + label + 1;
  printf("Read in eNodeB coordinate information done (%d).\n", *Tbs);
}

callingTimesCalculation(int UEInCoverage, FILE* parameter){
  int i, j;
  
 // int total = 0;
  for (i = 1; i <= UEInCoverage; i++){
    int nu = 1;
    int den = 1;
    j = 0;
    while (j < i){
      nu = nu * (UEInCoverage - j);
      den = den * (1 + j);
      j++;
    }
    if (nu > 0 && den > 0)
      fprintf(parameter, "combinations considering %d UEs: %d\n", i, (nu / den) * (int)pow(Modulation, i));
  }
  fprintf(parameter, "\n");
}
//------------------------------------------------------
//*******************************************************
void checkMSsites(int Tbs){
int i,j;
for(i=0;i<Ntp;i++){
	for(j=i+1;j<Ntp;j++)
		//if(i!=j)
			if(Xtp[i]==Xtp[j]&&Ytp[i]==Ytp[j]){
				printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		exit(1);
                }
                
  for(j=0;j<Tbs;j++)
    if(Xtp[i]==Xbs[j]&&Ytp[i]==Ybs[j]){
      printf("error!!MS[%d] and BS[%d] are at the same point.\n",i,j);
	exit(1);
    }  
}              
}
//------------------------------------------------------
//*******************************************************
void power(int Tbs, double TN){
	int i,j,k;

	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i<Tbs-mNbs)
					P[i][j][k]=pow(10,((-Gta-Gto-Gra-Gro+NF+FFM+BPL[j]+PL[i][j]+Sensitivity[j][k])/10))/1000;//+TN
				else
					P[i][j][k]=pow(10,((-mGta-mGto-Gra-Gro+NF+FFM+BPL[j]+PL[i][j]+Sensitivity[j][k])/10))/1000;
}
void PathLoss(double* ptrD, int* ptrPLT, int* ptrLT, int Tbs, double* ptrPL){
	int i,j;

	for(i=0;i<Tbs;i++){
		for(j=0;j<Ntp;j++){
			//printf("BS[%d]MS[%d] distance:%.2lf, location type:%d\n",i,j,D[i][j],LT[j]);
			if(LOSmin<=D[i][j]&&D[i][j]<NLOSmin){
				PathLossType[i][j]=1;	//LOS
				PL[i][j]=42.6+26*log10(D[i][j]/1000)+20*log10(BW/1000000);
				//printf("BS[%d]MS[%d] pathloss:%.2lf, type:%d\n",i,j,PL[i][j],PathLossType[i][j]);
			}
			else if(NLOSmin<=D[i][j]){
				PathLossType[i][j]=2;	//NLOS
				if(LT[j]==1)
					if(i<Tbs-mNbs)
						PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(BH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+((44.9-(6.55*log10(BH)))*(log10(D[i][j]/1000)))+3;
					else
						PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(mBH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+((44.9-(6.55*log10(mBH)))*(log10(D[i][j]/1000)))+3;
				else if(LT[j]==2||LT[j]==3)
					if(i<Tbs-mNbs)
						PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(BH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+(44.9-(6.55*log10(BH)))*(log10(D[i][j]/1000));
					else
						PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(mBH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+(44.9-(6.55*log10(mBH)))*(log10(D[i][j]/1000));
				else{
					printf("MS[%d] is in the unregular type.\n");
					exit(1);}
				//printf("BS[%d]MS[%d] pathloss:%.2lf, type:%d\n",i,j,PL[i][j],PathLossType[i][j]);
			}
			else{
				printf("BS[%d] and MS[%d] are at the same spot (%lf).\n",i,j, D[i][j]);
				exit(1);}
			//printf("BS[%d]MS[%d] DISTANCE:%.2lf, LT:%d, pathloss:%.2lf, type:%d\n",i,j,D[i][j],LT[j],PL[i][j],PathLossType[i][j]);
		}
	}
	return;
}
void FixedCellSize(int Tbs, int* DS, double* ptrD, double* P, double* ptrBR, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int* minM, int nUser){
	int i,j,k;
	double ShortestDistance=1000000000;
	int ConnectionAssignment[Ntp]={-1};
	int* ptrConnectionAssignment=&ConnectionAssignment[0];
	int mnbs=mNbs;
	double msbp=mSBP;
	double mabp=mABP;
	double sbp=SBP;
	double abp=ABP;
	
	for(j=0;j<Ntp;j++){
		for(i=0;i<Tbs;i++){
			if(D[i][j]<ShortestDistance){
				
				if(i<Tbs-mNbs&&D[i][j]<=Rbs){
					ShortestDistance=D[i][j];
					ConnectionAssignment[j]=i;
				}
				else if(i>=Tbs-mNbs&&D[i][j]<=rbs){
					ShortestDistance=D[i][j];
					ConnectionAssignment[j]=i;
				}
			}
		}
		if(ShortestDistance>Rbs){
			printf("%lf MS%d [%d][%d] locates outside the cells!\n Porgram Stops!!\n",ShortestDistance,j,Xtp[j],Ytp[j]);
			exit(1);
		}
		
		ShortestDistance=1000000000;
	}
	association(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ConnectionAssignment);
	zooming(Tbs, DS, minM, nUser, Modulation, Ntp, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, space, Xbs, Ybs, Xtp, Ytp, Rbs, rbs);
	switchOF(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ConnectionAssignment);
	//comparison4(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ConnectionAssignment);
	//comparison1(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ConnectionAssignment);
	int MConnectionAssignment[Ntp]={-1};
	int* ptrMConnectionAssignment=&MConnectionAssignment[0];
	for(j=0;j<Ntp;j++){
		for(i=0;i<Tbs-mNbs;i++){
			if(D[i][j]<ShortestDistance){
				if(D[i][j]<=Rbs){
					ShortestDistance=D[i][j];
					MConnectionAssignment[j]=i;
				}
			}
		}
		if(ShortestDistance>Rbs){
			printf("%lf MS%d [%d][%d] locates outside the cells!\n Porgram Stops!!\n",ShortestDistance,j,Xtp[j],Ytp[j]);
			exit(1);
		}
		ShortestDistance=1000000000;
	}
	int* MacroUserMod = malloc(Ntp*sizeof(int));
	//Macro(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrMConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, MacroUserMod);
	//free(ptrMConnectionAssignment);
	//printf("\n");
	//for(i=0;i<Ntp;i++)
	//  printf("%d\t",*(MacroUserMod+i));
	//NoMod(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, MacroUserMod);
	free(MacroUserMod);
	//BPSK(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ABP, SBP);
	//QAM64(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mP, minM, Rbs, rbs, space);
	//QAM16(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mP, minM, Rbs, rbs, space, ABP, SBP);
	//RANDOM(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mnbs, mP, minM, Rbs, rbs, space, ABP, SBP);
	//free(ptrConnectionAssignment);
	double TotalBandwidth=0.0, TotalPower=0.0, NOV=0.0;
	int TotalDS=0;
	double *Power=malloc(Tbs*sizeof(double));
	int *DataSub=malloc(Tbs*sizeof(int));
	for(i=0;i<Tbs;i++){
		*(Power+i)=0;
		*(DataSub+i)=0;
	}
	FILE *FixedCellSize;
 	if ((FixedCellSize=fopen("FixedCellSize.txt", "w")) == NULL){
		printf("\n\nerror!Fail to open file!");
		exit(1);
	}
	else
		printf("\n\nOpen FixedCellSize successfully!\n");

	fprintf(FixedCellSize,"\n#This file is meant to display the regular system.\n");
	fprintf(FixedCellSize,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	for(i=0;i<Tbs-mNbs;i++)
		fprintf(FixedCellSize," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],Rbs/space);
	fprintf(FixedCellSize,"\n\n");
	for(i=Tbs-mNbs;i<Tbs;i++)
		fprintf(FixedCellSize," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],rbs/space);
	fprintf(FixedCellSize,"\n\n");
	int UserModulation[Ntp]={0};
        for(i=0;i<Tbs;i++){
		fprintf(FixedCellSize,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
        	for(j=0;j<Ntp;j++)
		//	for(k=0;k<Modulation;k++)
				if(ConnectionAssignment[j]<0){
					printf("MS has no coonections! Program Stops!!\n");
					exit(1);
				}
				else if(ConnectionAssignment[j]==i){
				  
				  if (i < Tbs - mNbs){
					if(D[i][j]< (Rbs/Modulation))
						k=3;//k=6;
					else if(D[i][j]< (2 * (Rbs / Modulation)))
						k=2;//k=5;
					else if(D[i][j]< (3 * (Rbs / Modulation)))
						k=1;//k=4;
					//else if(D[i][j]<2800)
					//	k=2;//k=3;
					//else if(D[i][j]<3500)
						//k=1;//k=2;
					//else if(D[i][j]<4200)
					//	k=1;
					else
						k=0;
				  }
				  else{
				    	if(D[i][j]< (rbs/Modulation))
						k=3;//k=6;
					else if(D[i][j]< (2 * (rbs / Modulation)))
						k=2;//k=5;
					else if(D[i][j]< (3 * (rbs / Modulation)))
						k=1;//k=4;
					//else if(D[i][j]<2800)
					//	k=2;//k=3;
					//else if(D[i][j]<3500)
						//k=1;//k=2;
					//else if(D[i][j]<4200)
					//	k=1;
					else
						k=0;
				  }
//for(i=0;i<Ntp;i++){
					
		//for(j=0;j<Tbs;j++)
			//if(ConnectionAssignment[i]==j){
					Power[i]+=*(P+i*Ntp*Modulation+j*Modulation+k);
					if(i<Tbs-mNbs&&Power[i]>MP){
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					else if(i>=Tbs-mNbs&&Power[i]>mP){
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					DataSub[i]+=*(DS+j*Modulation+k);
						if(DataSub[i]>DSt){
							//printf("BS%d has too many subcarrier demands. Program Stops!\n",i);
							//exit(1);
						}
					
					TotalBandwidth+=BR[j];
					TotalDS+=*(DS+j*Modulation+k);
					//DataSub[j]+=*(DS+i);
					
					TotalPower+=*(P+i*Ntp*Modulation+j*Modulation+k);
			//}		
					NOV+=(*(BR+j)/(*(DS+j*Modulation+k)));
					
					fprintf(FixedCellSize,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
					fprintf(FixedCellSize,"BS[%d] %d %d MS[%d] %d %d DS %d power %lf mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),BR[j],k);
					UserModulation[j]=k;
	//}
				}
		//if(Power[i]>0){
			if(i<Tbs-mNbs){
				Power[i]+=ABP;
				TotalPower+=ABP;
			}
			else{
				Power[i]+=mABP;
				TotalPower+=mABP;
			}
		//}
// 		else{	
// 			if(i<Tbs-mNbs){	
// 				Power[i]+=SBP;
// 				TotalPower+=SBP;
// 			}
// 			else{
// 				Power[i]+=mSBP;
// 				TotalPower+=mSBP;
// 			}
// 		}
	}
	j=0;
	while(j<Modulation){
		fprintf(FixedCellSize,"\n\n");
		for(i=0;i<Ntp;i++)
			if(UserModulation[i]==j)
				fprintf(FixedCellSize,"MS[%d] %d %d DS %d BR %lf with modulation %d\n",i,Xtp[i],Ytp[i],*(DS+i*Modulation+j),BR[i],j);
		j++;
	}
	MCS(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, UserModulation);
	//comparison5(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, UserModulation);
	//comparison2(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR);
	fprintf(FixedCellSize,"\n\n");
	int temp;
	for(i=0;i<Tbs-mNbs;i++){
		temp=DataSub[i];
		for(j=Tbs-mNbs;j<Tbs;j++)
			if(*(minM+i*Tbs+j)==1)
				temp+=DataSub[j];
		if(temp>DSt){
			//printf("The group of BS%d has too many subcarrier demands. Program Stops!\n",i);
			//exit(1);
		}				
	}
	fprintf(FixedCellSize,"Total Power %lf\n",TotalPower);
        fprintf(FixedCellSize,"Total DSs: %d\n",TotalDS);
	fprintf(FixedCellSize,"Total BR: %lf\n",TotalBandwidth);
        int smallCellConnCount = 0;
	for(i=Tbs-mNbs;i<Tbs;i++){
	  for(j=0;j<Ntp;j++){
	    if(ConnectionAssignment[j]==i)
	      smallCellConnCount++;
	  }
	}
	fprintf(FixedCellSize,"Small Cell Connection Count: %d\n", smallCellConnCount);
	fprintf(FixedCellSize,"Spectrum-energy Efficiency: %lf\n",NOV/TotalPower);
	fprintf(FixedCellSize,"\n\n");
	
	fprintf(FixedCellSize," plot \"coordinates.txt\" index 0:0 using 2:3 title \"users\" with points, \"FixedCellSize.txt\" index 0:0 using 2:3:1 notitle with labels, \"FixedCellSize.txt\" index 0:0 using 2:3:4 title \"Macro Cell\" with circles, \"FixedCellSize.txt\" index 1:1 using 2:3:1 notitle with labels, \"FixedCellSize.txt\" index 1:1 using 2:3:4 title \"micro cell\" with circles");
	for(i=0;i<Tbs;i++)
		fprintf(FixedCellSize,", \"FixedCellSize.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i+2,i+2,i);
	for(i=0;i<Modulation;i++)
		fprintf(FixedCellSize,", \"FixedCellSize.txt\" index %d:%d using 2:3 title \"Modulation%d\" with points",Tbs+2+i,Tbs+2+i,i);
	fclose(FixedCellSize);
	OptModel(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment);
	EEMax(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment);
	SEMax(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment);
	SpecParti(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment);
	
	OptModelABS(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment);
	//printf("%lf", ABS1);
	//FixedABSpattern(1.0/8.0, 0.0, Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, maxCellSize, UserModulation, nUser, ptrAdjacentBSs);
	FixedABSpattern(1.0/8.0, 0.0, Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, maxCellSize, UserModulation, nUser, ptrAdjacentBSs);
	FixedABSpattern(2.0/8.0, 0.003, Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, maxCellSize, UserModulation, nUser, ptrAdjacentBSs);
	//FixedABSpattern(3.0/8.0, 0.003, Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, maxCellSize, UserModulation, nUser, ptrAdjacentBSs);
	//FixedABSpattern(3.0/8.0, 0.003, Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space, ptrD, ABP, SBP, maxCellSize, UserModulation, nUser, ptrAdjacentBSs);
	//OptFixedABSpattern(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mP, DSt, BR, ptrConnectionAssignment, 1.0/8.0, 0.0);
	binary_sleeping(Tbs, nUser, Ntp, Modulation, mNbs, DSt, mSBP, mABP, ABP, SBP, MP, mP, BP, 3.0 / 8.0);
	uniform_sleeping(Tbs, nUser, Ntp, Modulation, mNbs, DSt, mSBP, mABP, ABP, SBP, MP, mP, BP, 3.0 / 8.0);
	dynamic_activity_criterion_sleeping(Tbs, nUser, Ntp, Modulation, mNbs, DSt, mSBP, mABP, ABP, SBP, MP, mP, BP, 3.0 / 8.0);
	static_activity_criterion_sleeping(Tbs, nUser, Ntp, Modulation, mNbs, DSt, mSBP, mABP, ABP, SBP, MP, mP, BP, 3.0 / 8.0);
}
void DataSubcarriers(int* ptrDS, double* ptrBR, double* ptrDP){
	int j,k;

	for(j=0;j<Ntp;j++)
		for(k=0;k<Modulation;k++)
			  DS[j][k]=ceil(((1/N/3*9/8/4*5*FFT)/DP[j][k]*BR[j]/10)/12) * 1.5;
		  //DS[j][k]=ceil((1/N/3*9/8/4*5*FFT)/DP[j][k]*BR[j]/10);
	//printf("%.2lf\n",1/N/3*9/8*4/5*FFT);
		//BR[i]/(BW*1/N*DP*(1-1/5)*3/9)*8*FFT;
}
void LocationType(int* ptrLT, int* Xtp, int* Ytp){
	int j;

	for(j=0;j<Ntp;j++){
		if(X1urban<=Xtp[j]&&Xtp[j]<=X2urban&&Y1urban<=Ytp[j]&&Ytp[j]<=Y2urban)
			LT[j]=1;	//urban
		else if(X1suburban<=Xtp[j]&&Xtp[j]<=X2suburban&&Y1suburban<=Ytp[j]&&Ytp[j]<=Y2suburban)
			LT[j]=2;	//suburban
		else
			LT[j]=3;	//rural
	}
	return;
}
void bpl(int* ptrLT, int* ptrBPL){
	int i;
		
	for(i=0;i<Ntp;i++){
		if(LT[i]==1)
			BPL[i]=18;	//urban
		else if(LT[i]==2)
			BPL[i]=15;	//suburban
		else
			BPL[i]=12;	//rural
	}
	return;
}
void modulation(double TN, double* ptrSensitivity, double* ptrDP){//,double* ptrBR,  int* ptrModulation, double* ptrSNR){

	int j,k;

	for(j=0;j<Ntp;j++){
		DP[j][0]=1;
		DP[j][1]=2;
		DP[j][2]=3;
		//DP[j][3]=4.5;
		//DP[j][4]=4;
		//DP[j][5]=4.5;
		//DP[j][6]=4.5;

		Sensitivity[j][0]=-96;//TN+NF+9.4;
		Sensitivity[j][1]=-90;//TN+NF+11.2;
		Sensitivity[j][2]=-87;//TN+NF+16.4;
		//Sensitivity[j][3]=-80;//TN+NF+18.2;
		//Sensitivity[j][4]=TN+NF+22.7;
		//Sensitivity[j][5]=TN+NF+24.4;
		//Sensitivity[j][6]=TN+NF+24.4;
	}
// 	for(j=0;j<Ntp;j++){
// 		DP[j][0]=0.5;
// 		DP[j][1]=1;
// 		DP[j][2]=1.5;
// 		DP[j][3]=2;
// 		DP[j][4]=3;
// 		DP[j][5]=4;
// 		DP[j][6]=4.5;
// 
// 		Sensitivity[j][0]=TN+NF+6.4;
// 		Sensitivity[j][1]=TN+NF+9.4;
// 		Sensitivity[j][2]=TN+NF+11.2;
// 		Sensitivity[j][3]=TN+NF+16.4;
// 		Sensitivity[j][4]=TN+NF+18.2;
// 		Sensitivity[j][5]=TN+NF+22.7;
// 		Sensitivity[j][6]=TN+NF+24.4;
// 	}
	/*		
	for(i=1;i<=Ntp;i++){
		if(BR[i]<=0.86){
			SNR[i]=6.4;
			Modulation[i]=1;
			DP[i]=0.5;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=1.72){
			SNR[i]=9.4;
			Modulation[i]=2;
			DP[i]=1;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=2.58){
			SNR[i]=11.2;
			Modulation[i]=3;
			DP[i]=1.5;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=3.44){
			SNR[i]=16.4;
			Modulation[i]=4;
			DP[i]=2;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=5.16){
			SNR[i]=18.2;
			Modulation[i]=5;
			DP[i]=3;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=6.88){
			SNR[i]=22.7;
			Modulation[i]=6;
			DP[i]=4;
			Sensitivity[i]=TN+NF+SNR[i];}
		else {
			SNR[i]=24.4;
			Modulation[i]=7;
			DP[i]=4.5;
			Sensitivity[i]=TN+NF+SNR[i];}
	}*/
}
double ThermalNoise(){
	
return	-174+10*log10(BW*N*(720+120+1)/(double)1024);	//DSt+PS
	
//return tn;
}
void specifyTypes(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int i,j,k;

 	printf("binaries\n");
	fprintf(PowerSavingCPLEX,"binaries\n");
	fprintf(NoExistingUser,"binaries\n");
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				printf("a%d_%d_%d\n",i,j,k);//dBMi_j
 	for(i=0;i<Tbs;i++) 
		printf("b%d\n",i);		
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(NoExistingUser,"a%d_%d_%d\n",i,j,k);
 	for(i=0;i<Tbs;i++) 
		fprintf(NoExistingUser,"b%d\n",i);
	
	for(i=0;i<Tbs;i++)
  		for(j=0;j<Ntp;j++)
  			for(k=0;k<Modulation;k++)
   				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d\n",i,j,k);	
	for(i=0;i<Tbs;i++)
		fprintf(PowerSavingCPLEX,"dBA%d\n",i);
	
	fprintf(Feasibility, "binaries\n");
	for(i=0;i<Tbs;i++)
  		for(j=0;j<Ntp;j++)
  			for(k=0;k<Modulation;k++)
   				fprintf(Feasibility,"dBM%d_%d_%d\n",i,j,k);	
	for(i=0;i<Tbs;i++)
		fprintf(Feasibility,"dBA%d\n",i);
	fprintf(Feasibility, "end\n");
/*	for(i=1;i<=Tbs;i++) 
		fprintf(PowerSavingCPLEX,"dBS%d\n",i);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d\n",i,j);*/

return;
}
//*******************************************************
void bounds(int Tbs, FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int i, j,k;

	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				printf("0 <= u%d_%d_%d\n",i,j,k);//dBMi_j
 	for(i=0;i<Tbs;i++) 
		printf("0 <= v%d\n",i);
	//for(i=1;i<=Tbs;i++) 
		printf("0 <= t\n");		
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(NoExistingUser,"0 <= u%d_%d_%d\n",i,j,k);
 	for(i=0;i<Tbs;i++) 
		fprintf(NoExistingUser,"0 <= v%d\n",i);
	fprintf(NoExistingUser,"0 <= t\n");

	fprintf(PowerSavingCPLEX,"dBS=1\n");
	
	fprintf(Feasibility, "bounds\ndBS=1\n");
	//for(i=1;i<=Tbs;i++) 
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= dBM%d_%d\n",i,j);//dBMi_j
 	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBA%d\n",i);
	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBS%d <= 1\n",i);		
	//for(i=1;i<=Tbs;i++)
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= X%d_%d\n",i,j);
return;
}
   
   //*******************************************************
void checkAdjacentBSs(int* Xbs, int* Ybs, int* Xtp, int* Ytp, int Tbs) {
 // In this code, the candidate serving BSs for users are determined
 int i,j;

for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++)
		*(ptrAdjacentBSs+i*Ntp+j)=0;
 for(i=0;i<Tbs-mNbs;i++)
	for(j=0;j<Ntp;j++){
		
		if((space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])))<=MaxCellSize){
			//printf("Distance between user[%d] and BS[%d] is %lf\n",j,i,space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])));
			*(ptrAdjacentBSs+i*Ntp+j)=1;
		}
	}
for(i=Tbs-mNbs;i<Tbs;i++)
	for(j=0;j<Ntp;j++){
		
		if((space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])))<=maxCellSize){
			//printf("Distance between user[%d] and BS[%d] is %lf\n",j,i,space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])));
			*(ptrAdjacentBSs+i*Ntp+j)=1;
		}
	}/*
for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++)
		if(*(AdjacentBSs+i*Ntp+j)==1)
			printf("User[%d] is covered by BS[%d]\n",j,i);*/
}
//*******************************************************

void fillDistance(int* Xbs, int* Ybs, int* Xtp, int* Ytp, double* ptrD, int Tbs){
int i,j;
//double k=2;
for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++){
		D[i][j]=space*distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j]);
		//printf("%d %d %d %d D[%d][%d]=%.2lf\n",Xbs[i],Ybs[i],Xtp[j],Ytp[j],i,j,D[i][j]);
	}//printf("%.2lf",D[1][1]);
return;
}

void UECoordinates(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs, int nUser){
  int i, j, k, l;
  double ShortestDistance;
  double tempDistance = 0.0;
  srand(time(NULL));
  
/*  double UEacrossArea = 7.2;
double UEDensity = 425.69;
double gridSize = 10.0;
int UEhotSpotGrids = 4;*/
  double hotSpotUEDensity = UEDensity * 2;
  int hotSpotUENumber = hotSpotUEDensity * pow((gridSize * space / 1000), 2);
  int totalGrids = pow((map_y / space) / gridSize, 2);
  double nonHotSpotUEDensity = (Ntp - (UEhotSpotGrids * hotSpotUENumber)) / totalGrids;
  int axisGridNum = (map_y / space) / gridSize;
  for (i = 0; i < Ntp; i++){
  //for (i = totalGrids * nonHotSpotUEDensity - 1; i < Ntp; i++){
    if (i <= totalGrids * nonHotSpotUEDensity){
      *(Xtp + i) = rand() % (Px);
      *(Ytp + i) = rand() % (Py);
      *(ptrLT + i) = 3;
    }
    else if (i <= totalGrids * nonHotSpotUEDensity + 1 * hotSpotUENumber){
      *(Xtp + i) = rand() % (axisGridNum) + 62;
      *(Ytp + i) = rand() % (axisGridNum) + 86;
      //printf("UE%d is located in 1st hot spot.\n", i);
      *(ptrLT + i) = 2;
    }
    else if (i <= totalGrids * nonHotSpotUEDensity + 2 * hotSpotUENumber){
      *(Xtp + i) = rand() % (axisGridNum) + 118;
      *(Ytp + i) = rand() % (axisGridNum) + 61;
      //printf("UE%d is located in 2nd hot spot.\n", i);
      *(ptrLT + i) = 2;
    }
    else if (i <= totalGrids * nonHotSpotUEDensity + 3 * hotSpotUENumber){
      *(Xtp + i) = rand() % (axisGridNum) + 19;
      *(Ytp + i) = rand() % (axisGridNum) + 46;
      //printf("UE%d is located in 3rd hot spot.\n", i);
      *(ptrLT + i) = 2;
    }
    else if (i <= totalGrids * nonHotSpotUEDensity + 4 * hotSpotUENumber){
      *(Xtp + i) = rand() % (axisGridNum) + 98;
      *(Ytp + i) = rand() % (axisGridNum) + 21;
      //printf("UE%d is located in 4th hot spot.\n", i);
      *(ptrLT + i) = 2;
    }
    else{
      *(Xtp + i) = rand() % (Px);
      *(Ytp + i) = rand() % (Py);
      //printf("UE%d is located in 4th hot spot(%d %g %d).\n", i, totalGrids, nonHotSpotUEDensity, hotSpotUENumber);
      *(ptrLT + i) = 3;
    }
    ShortestDistance = 100000000000000;
    
    for (j = 0; j < Tbs - mNbs; j++){
      tempDistance = space * distance(Xbs[j], Ybs[j], Xtp[i], Ytp[i]);
      if (tempDistance < ShortestDistance)
	ShortestDistance = tempDistance;
    }
    if (ShortestDistance > Rbs)
      i = i - 1;
    for (k = i - 1; k > -1; k--)
      if (Xtp[i] == Xtp[k] && Ytp[i] == Ytp[k])
	i = i - 1;
      
    for (l = 0; l < Tbs; l++)
      if (Xbs[l] == Xtp[i] && Ybs[l] == Ytp [i])
	i = i - 1;
  }
  for(i=0;i<Tbs-mNbs;i++)
    for(j=0;j<Tbs;j++)
      *(ptrminM+i*Tbs+j)=0;
  *(ptrminM + 0 * Tbs + 29) = 1;
  *(ptrminM + 1 * Tbs + 29) = 1;
  *(ptrminM + 2 * Tbs + 29) = 1;
  *(ptrminM + 3 * Tbs + 29) = 1;
  *(ptrminM + 4 * Tbs + 28) = 1;
  *(ptrminM + 5 * Tbs + 28) = 1;
  *(ptrminM + 6 * Tbs + 30) = 1;
  *(ptrminM + 7 * Tbs + 31) = 1;
  *(ptrminM + 8 * Tbs + 30) = 1;
  *(ptrminM + 9 * Tbs + 30) = 1;
  *(ptrminM + 10 * Tbs + 32) = 1;
  *(ptrminM + 11 * Tbs + 31) = 1;
  *(ptrminM + 12 * Tbs + 34) = 1;
  *(ptrminM + 13 * Tbs + 32) = 1;
  //*(ptrminM + 13 * Tbs + 34) = 1;
  *(ptrminM + 14 * Tbs + 32) = 1;
  *(ptrminM + 15 * Tbs + 33) = 1;
  
  *(ptrminM + 16 * Tbs + 31) = 1;
  *(ptrminM + 17 * Tbs + 32) = 1;
 // *(ptrminM + 17 * Tbs + 34) = 1;
  *(ptrminM + 18 * Tbs + 31) = 1;
 // *(ptrminM + 18 * Tbs + 36) = 1;
  *(ptrminM + 19 * Tbs + 35) = 1;
  *(ptrminM + 20 * Tbs + 34) = 1;
  *(ptrminM + 21 * Tbs + 35) = 1;
  *(ptrminM + 22 * Tbs + 36) = 1;
  *(ptrminM + 23 * Tbs + 35) = 1;
  *(ptrminM + 24 * Tbs + 37) = 1;
  *(ptrminM + 25 * Tbs + 37) = 1;
  
  *(ptrminM + 26 * Tbs + 37) = 1;
  *(ptrminM + 27 * Tbs + 37) = 1;
  /* *(ptrminM + 28 * Tbs + 53) = 1;
  *(ptrminM + 28 * Tbs + 56) = 1;
  *(ptrminM + 29 * Tbs + 50) = 1;
  *(ptrminM + 29 * Tbs + 55) = 1;
  *(ptrminM + 30 * Tbs + 54) = 1;
  *(ptrminM + 31 * Tbs + 51) = 1;
  *(ptrminM + 31 * Tbs + 56) = 1;
  *(ptrminM + 32 * Tbs + 52) = 1;
  *(ptrminM + 32 * Tbs + 54) = 1;
  *(ptrminM + 33 * Tbs + 52) = 1;
  *(ptrminM + 33 * Tbs + 54) = 1;
  *(ptrminM + 34 * Tbs + 51) = 1;
  *(ptrminM + 34 * Tbs + 56) = 1;
  *(ptrminM + 35 * Tbs + 54) = 1;
  *(ptrminM + 36 * Tbs + 53) = 1;
  *(ptrminM + 36 * Tbs + 57) = 1;
  *(ptrminM + 37 * Tbs + 54) = 1;
  *(ptrminM + 38 * Tbs + 55) = 1;
  *(ptrminM + 39 * Tbs + 55) = 1;
  *(ptrminM + 40 * Tbs + 55) = 1;
  *(ptrminM + 41 * Tbs + 57) = 1;
  *(ptrminM + 42 * Tbs + 55) = 1;
  *(ptrminM + 43 * Tbs + 56) = 1;
  *(ptrminM + 44 * Tbs + 57) = 1;
  *(ptrminM + 45 * Tbs + 57) = 1;
  *(ptrminM + 46 * Tbs + 57) = 1;
  *(ptrminM + 47 * Tbs + 57) = 1;*/
  for (i = 0; i < Tbs; i++){
    //if (i >= Tbs - mNbs)
      *(ptrMm + i) = 0;
    if (i < Tbs - mNbs){
      for (j = 0; j < Tbs; j++){
	if (*(ptrminM + i * Tbs + j) == 1)
	  *(ptrMm + i) = *(ptrMm + i) + 1;
      }
      //printf("eNodeB%d has %d femto.\n", i, *(ptrMm + i));
    }
  }
  printf("UE Coordinate done.\n");
}

void fillCoordinatesTPs(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs, int nUser){
int i,j,k,l;

	double ShortestDistance;
	double tempDistance[Tbs];
	
	for(i=0;i<Tbs;i++)
		tempDistance[i]=0.0;
/*initialize random seed:*/
srand(time(NULL));

/*generate number:*/
for(i=0;i<nUser*UrbanUser;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2urban-X1urban+1)+X1urban;
	Ytp[i]=rand()%(Y2urban-Y1urban+1)+Y1urban;
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		
			i=i-1;
		
		//printf("%d renew\n",i);
		//continue;
	}
	//for(i=0;i<Ntp;i++)
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
        for (l=0; l<Tbs; l++)        
	  if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
	    i = i-1;
}
for(i=nUser*UrbanUser;i<nUser*(UrbanUser+SubUser);i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2suburban-X1suburban+1)+X1suburban;
	Ytp[i]=rand()%(Y2suburban-Y1suburban+1)+Y1suburban;
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		
			i=i-1;
		
		//printf("%d renew\n",i);
		//continue;
	}
	else if((Xtp[i]>X1urban&&Xtp[i]<X2urban)&&(Ytp[i]>Y1urban&&Ytp[i]<Y2urban)){
		i=i-1;
		//printf("%d renew\n",i);
	}
	
	else
	  for (l=0; l<Tbs; l++)        
	    if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
	      i = i-1;  	
	    
	for(k=i-1;k>-1;k--)
		//if(i!=j)
	  if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
	    i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
      
}
for(i=nUser*(UrbanUser+SubUser);i<nUser;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		
			i=i-1;
		
		//printf("%d renew\n",i);
		//continue;
	}
	else if((Xtp[i]>X1suburban&&Xtp[i]<X2suburban)&&(Ytp[i]>Y1suburban&&Ytp[i]<Y2suburban)){
		i=i-1;
		//printf("%d renew\n",i);
	}
	
	else
	  for (l=0; l<Tbs; l++)        
	    if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
	      i = i-1;
	    
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
                
}
for(i=nUser;i<nUser+ceil((Ntp-nUser)*UrbanUser);i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2urban-X1urban+1)+X1urban;
	Ytp[i]=rand()%(Y2urban-Y1urban+1)+Y1urban;
	
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		
			i=i-1;
		
		//printf("%d renew\n",i);
		//continue;
	}
	//for(i=0;i<Ntp;i++)
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
        
	for (l=0; l<Tbs; l++)        
	  if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
	    i = i-1;        
}

for(i=nUser+ceil((Ntp-nUser)*UrbanUser);i<nUser+ceil((Ntp-nUser)*(SubUser+UrbanUser));i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2suburban-X1suburban+1)+X1suburban;
	Ytp[i]=rand()%(Y2suburban-Y1suburban+1)+Y1suburban;
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		//printf("%d renew\n",i);
	}
	else if((Xtp[i]>X1urban&&Xtp[i]<X2urban)&&(Ytp[i]>Y1urban&&Ytp[i]<Y2urban)){
		i=i-1;
		//printf("%d renew\n",i);
	}
	else{
		for(k=i-1;k>-1;k--)
		//if(i!=j)
			if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
				i=i-1;
			
		for (l=0; l<Tbs; l++)        
		  if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
		    i = i-1;
	}
}
for(i=nUser+ceil((Ntp-nUser)*(SubUser+UrbanUser));i<Ntp;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		//printf("%d renew\n",i);
	}
	else if((Xtp[i]>X1suburban&&Xtp[i]<X2suburban)&&(Ytp[i]>Y1suburban&&Ytp[i]<Y2suburban)){
		i=i-1;
		//printf("%d renew\n",i);
	}
	else{
		for(k=i-1;k>-1;k--)
		//if(i!=j)
			if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
				i=i-1;
			
			
		for (l=0; l<Tbs; l++)        
		  if (Xbs[l]==Xtp[i] && Ybs[l]==Ytp[i])
		    i = i-1;
	}
}/*
for(i=41;i<=64;i++){
	Xtp[i]=(rand()%(376))+250;
	Ytp[i]=(rand()%(376))+250;
}*/
  printf("UE Coordinate done.\n");
return;
}

void fillCoordinatesBSs(int* Xbs, int* Ybs, int Gbs_x, int Gbs_y, int Nbs_x1, int Nbs_x2, int Nbs_y1, int Nbs_y2, int length, int Tbs, int* Mm, int* minM, int Px, int Py){
int Sx=1;
int Sy=1;
int Zy=Nbs_x1;
int Ay=Nbs_x1;
int By=Nbs_x2;
int Zx=Nbs_y1;
int Ax=Nbs_y1;
int Bx=Nbs_y2;
int i,j,l,k;
double m;
	double tDistance=0.0;
	int tBS=0;
	double sDistance=1000000000;
	srand(time(NULL));
 for(j=1; j<=(Gbs_y-1); j++) {
	 for(i=Sy; i<=Zy; i++)
		 Ybs[i-1] = j*length;
	 Sy=Zy+1;//printf("%d,%d\n",Sy,Zy);
	 if(j%2==1)
	 Zy+=By;
	 else
	 Zy+=Ay;	 
 }
 //printf("%d,%d",Xbs[4],Ybs[4]);
 //printf("%d",Tbs);
  //printf("Xrs[%d]=%lf, Yrs[%d]=%lf\n",i,Xrs[i],i,Yrs[i]);  
 for(m=1; m<=(Gbs_x-1); m+=1.5){
	for(l=Sx; l<=Sx+(Ay+By)*(Zx-1); l=l+Ay+By){
		Xbs[l-1]=m*length;
		//printf("BS[%d][%d]:(%d,%d)\n",l-1,l-1,Xbs[l-1],Ybs[l-1]);
	}
	if(((int)(m+2))%3==0){
		//printf("%lf %d\n",(int)(m+2),(int)(m+2));
		Sx+=Ay;
		Zx=Bx;
	}
	else{
		Sx=Sx-Ay+1;
		Zx=Ax;
	}
 }
	for(i=0;i<Tbs-mNbs;i++)
		for(j=0;j<Tbs;j++)
			*(minM+i*Tbs+j)=0;
	for(i=Tbs-mNbs;i<Tbs;i++){	//randomly generate
		sDistance=1000000000;
		Xbs[i]=rand()%(Px+1);
		Ybs[i]=rand()%(Py+1);
		//do{
		//  Xbs[i]=rand()%(Px+1);
		//  Ybs[i]=rand()%(Py+1);
		//}while(Xbs[i]<X1urban || Xbs[i]>X2urban || Ybs[i]<Y1urban || Ybs[i]>Y2urban);
		for(j=0;j<Tbs-mNbs;j++){
			tDistance=space*distance(Xbs[j],Ybs[j],Xbs[i],Ybs[i]);
			if(tDistance<sDistance){
				sDistance=tDistance;
				tBS=j;
			}
		}
		l=0;
		if (sDistance>Rbs){
		  i=i-1;
		  l=1;
		  printf("renew BS location\n");
		}
		for(k=i-1;k>-1;k--)
		//if(i!=j)
			if(Xbs[i]==Xbs[k]&&Ybs[i]==Ybs[k]){
				i=i-1;
				l=1;
			}
			
		if(l==0){
			*(Mm+tBS)+=1;
			*(minM+tBS*Tbs+i)=1;
			//printf("BS%d has %d micro BSs.\n",tBS,*(Mm+tBS));
		}	
		
		
	}
	j=0;
	for(i=0;i<Tbs;i++){
		j+=*(Mm+i);
		printf("BS%d has %d micro BSs.\n",i,*(Mm+i));
	}
	if(j!=mNbs){
		printf("Micro Cell Number error! %d\n",j);
		exit(1);
	}
 return;
}
//*******************************************************

double distance(double xa, double ya, double xb, double yb) {
 double value;
 //double temp1, temp2;

// temp1 = xa-xb;
// temp1 = pow(temp1,2);
// temp2 = ya-yb;
// temp2 = pow(temp2,2);
// value = sqrt(temp1+temp2);

 value = sqrt( pow((xa-xb),2) + pow((ya-yb),2) );

 return value;
}
//*******************************************************
fillMs(double* ptrBR, int* Xtp, int* Ytp) { // Mbr,Mbt,Mrr,Mrt/* The grid is a square with length=width. The BS is at the top line in the middle.
/*length is even, the BS will coincide with a RS site. Then we "cancel" this RS,
which means mbr=0. See function_fillMbr.ps illustration.*/
	int i;
	//double dist;
	double j;
	srand(time(NULL));

	for(i=0;i<Ntp;i++){
		//BR[i]
		j=(rand()%(MBR+1));
		while(j==0)
			j=(rand()%(MBR+1));
		BR[i]=j/100;
		//printf("%.2lf\n",j);
		//printf("%.2lf\n",BR[i]);
		if (61 <= *(Xtp + i) <= 77 && 85 <= *(Ytp + i) <= 101){
		  BR[i] = BR[i] * 1.5;
		  *(ptrLT + i) = 1;
		}
	}
	return;
}
void constraint11(int Tbs,FILE* PowerSavingCPLEX,int nUser, FILE* Feasibility){
	int i,j,k;

	
	for(j=nUser;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(k==Modulation-1&&i==Tbs-1)
					printf("u%d_%d_%d - t = 0\n",i,j,k);
				else
					printf("u%d_%d_%d + ",i,j,k);
	for(j=nUser;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d = 1\n",i,j,k);
				else
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d + ",i,j,k);
				
	for(j=nUser;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&k==Modulation-1)
					fprintf(Feasibility,"dBM%d_%d_%d = 1\n",i,j,k);
				else
					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
}
void constraint10(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility){
	int i,j,k;

	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(k==Modulation-1)
					printf("u%d_%d_%d ",i,j,k);
				else
					printf("u%d_%d_%d + ",i,j,k);
			}
			printf("- %d t <= 0\n",AdjacentBSs[i][j]);
		}
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(k==Modulation-1)
					fprintf(NoExistingUser,"u%d_%d_%d ",i,j,k);
				else
					fprintf(NoExistingUser,"u%d_%d_%d + ",i,j,k);
			}
			fprintf(NoExistingUser,"- %d t <= 0\n",AdjacentBSs[i][j]);
		}
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(k==Modulation-1)
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d ",i,j,k);
				else
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d + ",i,j,k);
			}
			fprintf(PowerSavingCPLEX,"<= %d\n",AdjacentBSs[i][j]);
		}
		
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(k==Modulation-1)
					fprintf(Feasibility,"dBM%d_%d_%d ",i,j,k);
				else
					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
			}
			fprintf(Feasibility,"<= %d\n",AdjacentBSs[i][j]);
		}
}
//*******************************************************
void constraint9(int Tbs, FILE* NoExistingUser, FILE* Feasibility) {
 int i,j,k;

	for(i=0;i<Tbs-mNbs;i++){
		printf("%lf v%d + ",ABP-SBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d",P[i][j][k],i,j,k);
				else
					printf("%lf u%d_%d_%d + ",P[i][j][k],i,j,k);
			
	}
	for(i=Tbs-mNbs;i<Tbs;i++){
		printf(" + %lf v%d",mABP-mSBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
					printf(" + %lf u%d_%d_%d",P[i][j][k],i,j,k);
				else
					printf(" + %lf u%d_%d_%d",P[i][j][k],i,j,k);
			
	}
	printf(" + %lf t = 1\n",((Tbs-mNbs)*SBP)+(mNbs*mSBP));
	for(i=0;i<Tbs-mNbs;i++){
		fprintf(NoExistingUser,"%lf v%d + ",ABP-SBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(NoExistingUser,"%lf u%d_%d_%d",P[i][j][k],i,j,k);
				else
					fprintf(NoExistingUser,"%lf u%d_%d_%d + ",P[i][j][k],i,j,k);
			
	}
	for(i=Tbs-mNbs;i<Tbs;i++){
		fprintf(NoExistingUser," + %lf v%d",mABP-mSBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(NoExistingUser," + %lf u%d_%d_%d",P[i][j][k],i,j,k);
				else
					fprintf(NoExistingUser," + %lf u%d_%d_%d",P[i][j][k],i,j,k);
			
	}
	fprintf(NoExistingUser," + %lf t = 1\n",((Tbs-mNbs)*SBP)+(mNbs*mSBP));

	for(j=0;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(Feasibility,"%lf dBM%d_%d_%d",BR[j]/DS[j][k],i,j,k);
				else 
					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",BR[j]/DS[j][k],i,j,k);
	fprintf(Feasibility," - ");
	for(i=0;i<Tbs-mNbs;i++){
		fprintf(Feasibility,"%lf dBA%d - ",(ABP-SBP)*limit,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(Feasibility,"%lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
				else
					fprintf(Feasibility,"%lf dBM%d_%d_%d - ",P[i][j][k]*limit,i,j,k);
	}
	for(i=Tbs-mNbs;i<Tbs;i++){
		fprintf(Feasibility," - %lf dBA%d",(mABP-mSBP)*limit,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(Feasibility," - %lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
				else
					fprintf(Feasibility," - %lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
			
	}
	fprintf(Feasibility," - %lf dBS >= 0\n",(((Tbs-mNbs)*SBP)+(mNbs*mSBP))*limit);

 return;
}
//*******************************************************
void constraint8(int Tbs, FILE* NoExistingUser) {
 int i,j,k;
	int Q2=1002;

	for(i=0;i<Tbs;i++){
		printf("v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		printf("v%d - t + %d b%d <= %d\n",i,Q2,i,Q2);
		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d <= %d\n",i,Q2,i,Q2);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				printf("u%d_%d_%d - t - %d a%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
				printf("u%d_%d_%d - t + %d a%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
				//printf("w%d_%d_%d - t - %d c%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
				//printf("w%d_%d_%d - t - %d c%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
				//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
				//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d <= %d\n",i,j,Q2,i,j,Q2);
				//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
				//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d <= %d\n",i,j,Q2,i,j,Q2);
			}
	}
	for(i=0;i<Tbs;i++){
		fprintf(NoExistingUser,"v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		fprintf(NoExistingUser,"v%d - t + %d b%d <= %d\n",i,Q2,i,Q2);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				fprintf(NoExistingUser,"u%d_%d_%d - t - %d a%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
				fprintf(NoExistingUser,"u%d_%d_%d - t + %d a%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
			}
	}
 return;
}
//*******************************************************
void constraint7(int Tbs, FILE* NoExistingUser) {
	int i,j,k;
	int Q1=1002;

	for(i=0;i<Tbs;i++){
		fprintf(NoExistingUser,"v%d - %d b%d <= 0\n",i,Q1,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				fprintf(NoExistingUser,"u%d_%d_%d - %d a%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
			}
	}
	for(i=0;i<Tbs;i++){
		printf("v%d - %d b%d <= 0\n",i,Q1,i);
		//fprintf(PowerSavingCPLEX,"v%d <= %d t\n",i,Q1);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				printf("u%d_%d_%d - %d a%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
				//printf("w%d_%d_%d - %d c%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
				//fprintf(PowerSavingCPLEX,"u%d_%d <= %d t\n",i,j,Q1);
				//fprintf(PowerSavingCPLEX,"w%d_%d <= %d t\n",i,j,Q1);
			}
	}
return;
}
//*******************************************************
void constraint6(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int i,j,k;
	int Q=1002;

	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				//if(j==Ntp-1&&k==Modulation-1)
					printf("u%d_%d_%d - v%d <= 0\n",i,j,k,i);
				/*else
					printf("u%d_%d_%d + ",i,j,k);*/
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				//if(j==Ntp-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d - dBA%d <= 0\n",i,j,k,i);
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(NoExistingUser,"u%d_%d_%d - v%d <= 0\n",i,j,k,i);
			
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(Feasibility,"dBM%d_%d_%d - dBA%d <= 0\n",i,j,k,i);
 return;
}
//*******************************************************
void constraint5(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int i,j,k;

 	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					printf("u%d_%d_%d - v%d >= 0\n",i,j,k,i);
				else
					printf("u%d_%d_%d + ",i,j,k);
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(NoExistingUser,"u%d_%d_%d - v%d >= 0\n",i,j,k,i);
				else
					fprintf(NoExistingUser,"u%d_%d_%d + ",i,j,k);
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d - dBA%d >= 0\n",i,j,k,i);
				else
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d + ",i,j,k);
				
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(Feasibility,"dBM%d_%d_%d - dBA%d >= 0\n",i,j,k,i);
				else
					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
return;
}
//*******************************************************
void constraint4(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
 int i,j,k;

 for(i=0; i<Tbs; i++) 
  //printf("dbt%d",i);
	 for(j=0; j<Ntp; j++)
		for(k=0;k<Modulation;k++)
			if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
				printf("u%d_%d_%d",i,j,k);
			else 
				printf("u%d_%d_%d +",i,j,k);
	 
  
 
printf(" - %.2lf t >= 0\n",Ntp*(1-BP));
	for(i=0; i<Tbs; i++) 
	 	for(j=0; j<Ntp; j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(NoExistingUser,"u%d_%d_%d",i,j,k);
				else 
					fprintf(NoExistingUser,"u%d_%d_%d +",i,j,k);
	fprintf(NoExistingUser," - %.2lf t >= 0\n",Ntp*(1-BP));
	for(i=0; i<Tbs; i++) 
  //printf("dbt%d",i);
	 for(j=0; j<Ntp; j++)
		for(k=0;k<Modulation;k++)
			if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d",i,j,k);
			else 
				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d +",i,j,k);
	 
  
 
fprintf(PowerSavingCPLEX," >= %.2lf\n",Ntp*(1-BP));

	for(i=0; i<Tbs; i++) 
		for(j=0; j<Ntp; j++)
			for(k=0; k<Modulation; k++)
				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(Feasibility,"dBM%d_%d_%d",i,j,k);
				else 
					fprintf(Feasibility,"dBM%d_%d_%d +",i,j,k);
	fprintf(Feasibility," >= %.2lf\n",Ntp*(1-BP));
 return;
}
//*******************************************************
void constraint3(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, int nUser, FILE* Feasibility) {
 int i,j,k;
	 for(i=0; i<nUser; i++)
		 for(j=0; j<Tbs; j++) 
		 	for(k=0;k<Modulation;k++)
				if(j==Tbs-1&&k==Modulation-1)
					printf("u%d_%d_%d - t <= 0\n",j,i,k);
				else
					printf("u%d_%d_%d +",j,i,k);
	for(i=0; i<Ntp; i++)
		 for(j=0; j<Tbs; j++) 
		 	for(k=0;k<Modulation;k++)
				if(j==Tbs-1&&k==Modulation-1)
					fprintf(NoExistingUser,"u%d_%d_%d - t <= 0\n",j,i,k);
				else
					fprintf(NoExistingUser,"u%d_%d_%d +",j,i,k);
	for(i=0; i<nUser; i++)
		 for(j=0; j<Tbs; j++) 
		 	for(k=0;k<Modulation;k++)
				if(j==Tbs-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d <= 1\n",j,i,k);
				else
					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d +",j,i,k);
	
	for(i=0; i<nUser; i++)
		 for(j=0; j<Tbs; j++) 
		 	for(k=0;k<Modulation;k++)
				if(j==Tbs-1&&k==Modulation-1)
					fprintf(Feasibility,"dBM%d_%d_%d <= 1\n",j,i,k);
				else
					fprintf(Feasibility,"dBM%d_%d_%d +",j,i,k);
 return;
}

//*******************************************************
void constraint2(int Tbs,FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int i,j,k;

	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d", P[i][j][k],i, j,k);
				else 
					printf("%lf u%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	printf(" - %lf t <= 0\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d", P[i][j][k],i, j,k);
				else 
					printf("%lf u%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	printf(" - %lf t <= 0\n",mP);
	}
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(NoExistingUser,"%lf u%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(NoExistingUser,"%lf u%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	fprintf(NoExistingUser," - %lf t <= 0\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(NoExistingUser,"%lf u%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(NoExistingUser,"%lf u%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	fprintf(NoExistingUser," - %lf t <= 0\n",mP);
	}
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	fprintf(PowerSavingCPLEX," <= %lf\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
		}
	fprintf(PowerSavingCPLEX," <= %lf\n",mP);
	}
	
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(Feasibility,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
		}
		fprintf(Feasibility," <= %lf\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(Feasibility,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
				else 
					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
		}
		fprintf(Feasibility," <= %lf\n",mP);
	}
 return;
}//*******************************************************

void constraint1(int Tbs, int* ptrDS, FILE* PowerSavingCPLEX, int* Mm, int* minM, FILE* NoExistingUser, FILE* Feasibility) {
	int i,j,k,l,m,n;
	//int o=0;
	for(i=0; i<Tbs-mNbs; i++){
		//o+=*(Mm+i);
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					printf("%d u%d_%d_%d", DS[j][k],i, j,k);
					//if(*(Mm+i)>0)
					for(l=Tbs-mNbs;l<Tbs;l++)//for(l=Tbs-mNbs+o-*(Mm+i);l<Tbs-mNbs+o;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									//if(m==Ntp-1&&n==Modulation-1)
										printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
									//else
										//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
				}		
				else 
					printf("%d u%d_%d_%d + ",DS[j][k], i, j,k);
		}
	printf(" - %d t <= 0\n",DSt);	//DSt*200 for 1 second
	}
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(NoExistingUser,"%d u%d_%d_%d", DS[j][k],i, j,k);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(NoExistingUser," + %d u%d_%d_%d",DS[m][n],l,m,n);
				}		
				else 
					fprintf(NoExistingUser,"%d u%d_%d_%d + ",DS[j][k], i, j,k);
		}
	fprintf(NoExistingUser," - %d t <= 0\n",DSt);	
	}
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(PowerSavingCPLEX," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
				}		
				else 
					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
		}
	fprintf(PowerSavingCPLEX," <= %d\n",DSt);	
	}
	
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(Feasibility,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(Feasibility," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
				}		
				else 
					fprintf(Feasibility,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
		}
	fprintf(Feasibility," <= %d\n",DSt);	
	}	
 return;
}//*******************************************************

void objective(int Tbs, FILE* PowerSavingCPLEX, FILE* NoExistingUser, FILE* Feasibility) {
	int count,i,k;
 
	printf("maximize\n");
	fprintf(PowerSavingCPLEX,"minimize\n");
	fprintf(NoExistingUser,"maximize\n");
	for(count=0;count<Ntp;count++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d",BR[count]/DS[count][k],i,count,k);
				else 
					printf("%lf u%d_%d_%d + ",BR[count]/DS[count][k],i,count,k);
	for(count=0;count<Ntp;count++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(NoExistingUser,"%lf u%d_%d_%d",BR[count]/DS[count][k],i,count,k);
				else 
					fprintf(NoExistingUser,"%lf u%d_%d_%d + ",BR[count]/DS[count][k],i,count,k);
	for(count=0; count<Tbs-mNbs; count++){
		for(i=0;i<Ntp;i++)
			for(k=0;k<Modulation;k++)
				if(count==0&&i==0&&k==0)
					fprintf(PowerSavingCPLEX,"%lf dBS + %lf dBM%d_%d_%d + ",SBP*(Tbs-mNbs),P[count][i][k],count,i,k);
				else if(count!=Tbs-mNbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",ABP-SBP,count);
					//fprintf(PowerSavingCPLEX," + %lf + ",SBP);
				}
				else if(count==Tbs-mNbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",ABP-SBP,count);						
											
					
					//fprintf(PowerSavingCPLEX," + %lf\n",SBP*Tbs);
				} 
				else
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[count][i][k],count,i,k);
		
	}		
	for(count=Tbs-mNbs; count<Tbs; count++)
		for(i=0;i<Ntp;i++)
			for(k=0;k<Modulation;k++)
				if(count==Tbs-mNbs&&i==0&&k==0)
					fprintf(PowerSavingCPLEX,"%lf dBS + %lf dBM%d_%d_%d + ",mSBP*mNbs,P[count][i][k],count,i,k);
				else if(count!=Tbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",mABP-mSBP,count);
					//fprintf(PowerSavingCPLEX," + %lf + ",SBP);
				}
				else if(count==Tbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d",mABP-mSBP,count);
					//fprintf(PowerSavingCPLEX," + %lf\n",SBP*Tbs);
				} 
				else
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[count][i][k],count,i,k);
				
	fprintf(Feasibility, "maximize\n0\nst\n");
return;
}
