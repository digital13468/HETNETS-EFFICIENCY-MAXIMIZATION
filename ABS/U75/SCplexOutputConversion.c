#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"header/CoverageComputation.h"
#include<math.h>

#define Tbs 38
#define Ntp 1368
#define BP 0.075
#define SBP 75.00	//Power for BS in sleep mode
#define ABP 130.0	//Basic power for BS in active mode
#define Modulation 3
#define mNbs 10
#define mSBP 2.8	//Power for micro BSs in sleep mode
#define mABP 4.8	//Basic power for micro BSs in active mode
int Xb[Tbs]={0};
int Yb[Tbs]={0};
int B[Tbs]={0};
int Xtp[Ntp]={0};
int Ytp[Ntp]={0};
int dBM[Tbs][Ntp][Modulation];
int *ptrdBM;
int dBA[Tbs]={0};
int A[Tbs][Ntp][Modulation];
int C[Tbs][Ntp][Modulation];
double P[Tbs][Ntp][Modulation];
double BR[Ntp]={0};
double u[Tbs][Ntp][Modulation];
double w[Tbs][Ntp][Modulation];
int DataSub[Ntp][Modulation];
double CumulativeP[Tbs]={0};
double CumulativeBR[Tbs]={0};
double v[Tbs]={0};
int CumulativeDS[Tbs]={0};
int CumSerMS[Tbs]={0};
double TotalP=0.0;
double TotalBR=0.0;
double t=0.0;
double NOV=0.0;
int TotalDS=0;
int TotSerMS=0;

double x[Tbs]={0};

/*InfoPrinter(FILE* CplexOut){
               int i,j;
               for(i=0;i<Tbs;i++)
               fprintf(CplexOut,"       BS%3d",i+1);
               fprintf(CplexOut,"\n");
               
               for(i=0;i<Ntp;i++){
                                  fprintf(CplexOut,"MS%4d%4d",i+1,dBM[0][i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"      %6d",dBM[j][i]);
                                  fprintf(CplexOut,"\n");
                                  
                                  fprintf(CplexOut," Power%10.5f",P[0][i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"  %10.5f",P[j][i]);
                                  fprintf(CplexOut,"\n");
                                  
                                  fprintf(CplexOut," DS %6d",DataSub[i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"      %6d",DataSub[i]);
                                  fprintf(CplexOut,"\n");
                                  }
}*/

double * ConnectionPrinter(FILE* CplexOut){
               int i,j,k;
	double *coverage;
	
	coverage=CoverageComputation(ptrdBM,Tbs,Ntp,Xb,Yb,Xtp,Ytp,Modulation);

fprintf(CplexOut,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	fprintf(CplexOut,"#plot \"coordinates.txt\" index 0:0 using 2:3 title \"users\" with points, \"SpectrumEfficiencyCPLEX_Output.txt\" index 0:0 using 2:3:1 notitle with labels, \"SpectrumEfficiencyCPLEX_Output.txt\" index 0:0 using 2:3:4 title \"Macro Cell\" with circles, \"SpectrumEfficiencyCPLEX_Output.txt\" index 1:1 using 2:3:1 notitle with labels, \"SpectrumEfficiencyCPLEX_Output.txt\" index 1:1 using 2:3:4 title \"micro cell\" with circles");
	for(i=0;i<Tbs;i++)
		fprintf(CplexOut,", \"SpectrumEfficiencyCPLEX_Output.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i+2,i+2,i);
	for(i=0;i<Modulation;i++)
		fprintf(CplexOut,", \"SpectrumEfficiencyCPLEX_Output.txt\" index %d:%d using 2:3 title \"Modulation%d\" with points",Tbs+2+i,Tbs+2+i,i);
// fprintf(CplexOut,"# plot \"SpectrumEfficiencyCPLEX_Output.txt\" index 0:0 using 2:3:1 notitle with labels, \"SpectrumEfficiencyCPLEX_Output.txt\" index 1:1 using 5:6:(($13)*0.001) title \"Spectrum-energy efficiency\" with circles fill solid lt rgb \"cyan\"");// \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
/*for(i=1;i<=Tbs;i++)
fprintf(CplexOut,", \"outfile.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i,i,i);*/
fprintf(CplexOut,"\n\n");

for(i=0;i<Tbs-mNbs;i++)
	fprintf(CplexOut," \"BS[%d]\" %d %d %g\n",i,Xb[i],Yb[i],*(coverage+i));
fprintf(CplexOut,"\n\n");
for(i=Tbs-mNbs;i<Tbs;i++)
	fprintf(CplexOut," \"BS[%d]\" %d %d %g\n",i,Xb[i],Yb[i],*(coverage+i));
fprintf(CplexOut,"\n\n");
for(i=0;i<Tbs;i++){
	fprintf(CplexOut,"BS[%d] %d %d BS[%d] %d %d\n",i,Xb[i],Yb[i],i,Xb[i],Yb[i]);
        for(j=0;j<Ntp;j++)
        	for(k=0;k<Modulation;k++)
			if(dBM[i][j][k]==1){
				//fprintf(CplexOut,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xb[i],Yb[i],i+1,Xb[i],Yb[i]);
				fprintf(CplexOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW SE %g with modulation %d\n",i,Xb[i],Yb[i],j,Xtp[j],Ytp[j],DataSub[j][k],P[i][j][k],BR[j]/DataSub[j][k]/P[i][j][k],k);
			}

	fprintf(CplexOut,"\n\n");
}
	j=0;
	while(j<Modulation){
		fprintf(CplexOut,"\n\n");
		for(i=0;i<Tbs;i++)
			for(k=0;k<Ntp;k++)
				if(dBM[i][k][j]==1)
					fprintf(CplexOut,"MS[%d] %d %d DS %d power %g mW select Mod %d\n",k,Xtp[k],Ytp[k],DataSub[k][j],P[i][k][j],j);
		j++;
	}	
	fprintf(CplexOut,"\n\n");
}
                                  
main(){
       
FILE *infile, *outfile; 
int i,j,k;
char id[1000];
char tempC;
int a, b, c;
int BS,MS,value,Mod;
double fvalue;
size_t loc;
char string[12]="<variables>" ;
char string1[13]="</variables>" ;
char string2[10]="<variable";
ptrdBM=&dBM[0][0][0];
if ((infile=fopen("SECRD.txt", "r")) == NULL)
	printf("Fail to open file!");
else
	printf("Open file successfully!\n");
/*if((outfile=fopen("outfile.txt","w"))==NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");*/
	for(j=0;j<Ntp;j++)
		for(k=0;k<Modulation;k++){
			DataSub[j][k]=0;
			for(i=0;i<Tbs;i++){
				dBM[i][j][k]=0;
				P[i][j][k]=0;
				A[i][j][k]=0;
				C[i][j][k]=0;
				u[i][j][k]=0;
				w[i][j][k]=0;
			}
		}
	while (fscanf(infile, "%s", &id)!=EOF){
		if (strcmp(id,string)==0){
			fscanf(infile, "%s", &id);
			while(strcmp(id,string1)!=0){//printf("Outer Scanning\n");
				while(strcmp(id,string2)==0){
					fscanf(infile, " name=\"");
					fscanf(infile, "%c", &tempC);
					//printf("%c\n",tempC);
					if(tempC=='a'){
						fscanf(infile, "%d_%d_%d\"", &BS, &MS, &Mod);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						if(fvalue<=0.5)
						  A[BS][MS][Mod]=0;
						else
						  A[BS][MS][Mod]=1;
						//fscanf(infile, "%s",&id);
						//printf("%s\n",id);
						//fscanf(infile, "%s",&id);
						//printf("%s\n",id);
						//fscanf(infile, "%s",&id);
						//printf("%s\n",id);
					}
					else if(tempC=='b'){
						fscanf(infile, "%d\"", &BS);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						if(fvalue<=0.5)
						  B[BS]=0;
						else
						  B[BS]=1;
					}
					else if(tempC=='c'){
						fscanf(infile, "%d_%d_%d\"", &BS, &MS, &Mod);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%d\"/>", &value);
						C[BS][MS][Mod]=value;
						printf("%s\n",id);
					}
					else if(tempC=='u'){
						fscanf(infile, "%d_%d_%d\"", &BS, &MS, &Mod);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						//printf("%lf\n",fvalue);
						u[BS][MS][Mod]=fvalue;
						//printf("%lf\n",fvalue);
					}
					else if(tempC=='v'){
						fscanf(infile, "%d\"", &BS);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						v[BS]=fvalue;
					}
					else if(tempC=='w'){
						fscanf(infile, "%d_%d_%d\"", &BS, &MS, &Mod);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						w[BS][MS][Mod]=fvalue;
					}
					else if(tempC=='t'){
						fscanf(infile, "\"");
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						t=fvalue;
					}
					
					else if(tempC=='x'){
						fscanf(infile, "%d\"", &BS);
						fscanf(infile, "%s", &id);
						fscanf(infile, " value=\"%lf\"/>", &fvalue);
						x[BS]=fvalue;
					}
					/*else if(tempC=='d'){
						fscanf(infile, "%c%c", &tempC,&tempC);
						//printf("%c\n",tempC);
						if(tempC=='M'){
							fscanf(infile, "%d_%d\"", &BS, &MS);
							fscanf(infile, "%s", &id);
							fscanf(infile, " value=\"%d\"/>", &value);
							dBM[BS-1][MS-1]=value;
							printf("%d %d\n",BS,MS);
						}
						else{
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
						}
					}*/
					//printf("%s\n",id);
					fscanf(infile, "%s", &id);
				/*fscanf(infile, " name=\"dBM%d_%d\"", &BS, &MS);
				fscanf(infile, "%s", &id);
				fscanf(infile, " value=\"%d\"/>", &value);
				dBM[BS-1][MS-1]=value;
				fscanf(infile, "%s", &id);
					//printf("Inner Scanning\n");
					if(BS==Tbs&&MS==Ntp)
						strcpy(id,string1);*/
				}                                                          
			}                                                            
		}
	
	}

/*for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
fprintf(outfile,"dBM%d_%d=%d\n",i+1,j+1,dBM[i][j]);*/
fclose(infile);
//fclose(outfile);

FILE *CplexOut;
FILE *CplexIn;
//int BS,MS,value;
int Xbs, Ybs, Xms, Yms, DS;
double p,br;
char id1[1000];

if ((CplexIn=fopen("outfile1.txt", "r")) == NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");
if ((CplexOut=fopen("SpectrumEfficiencyCPLEX_Output.txt", "w")) == NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");
fscanf(CplexIn, "%s", &id1);
while (fscanf(CplexIn, "%s", &id1)!=EOF){
      //fprintf(CplexOut, "%s\n",id1);
      fscanf(CplexIn, " BS[%d]=%d,%d", &BS, &Xbs, &Ybs);
	Xb[BS]=Xbs;
	Yb[BS]=Ybs;
      //fprintf(CplexOut,"BS[%d]=(%d,%d) ", BS, Xbs, Ybs);
      fscanf(CplexIn, " MS[%d]=%d,%d", &MS, &Xms, &Yms);
	Xtp[MS]=Xms;
	Ytp[MS]=Yms;
      //fprintf(CplexOut,"MS[%d]=(%d,%d) ", MS, Xms, Yms);
      fscanf(CplexIn, " with modulation %d", &Mod);
      fscanf(CplexIn, " DS=%d", &DS);
      DataSub[MS][Mod]=DS;
      //fprintf(CplexOut,"DS:%d ", DS);
      fscanf(CplexIn, " power=%lf mW", &p);
      P[BS][MS][Mod]=p;
	fscanf(CplexIn," BW=%lf",&br);
	BR[MS]=br;
      //fprintf(CplexOut,"power:%g\n", p);
      fscanf(CplexIn," %s",&id1);

}
/*for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
fprintf(outfile,"dBM%d_%d=%d P:%g DS:%d\n",i+1,j+1,dBM[i][j],P[i][j],DataSub[j]);*/
fclose(CplexIn);

for(i=0;i<Tbs;i++){
  if(v[i]/t >= 0.5)
	dBA[i]=1;
  else
    dBA[i]=0;
	for(j=0;j<Ntp;j++){
		for(k=0;k<Modulation;k++){
		//fprintf(CplexOut,"%d %d %lf %lf %lf\n",i,j,u[i][j], t, u[i][j]/t);
		  if(u[i][j][k]/t>=0.5)
			dBM[i][j][k]=1;
		  else
		    dBM[i][j][k]=0;
		}
	}
}
		
for(i=0;i<Tbs;i++){
	if(dBA[i]!=B[i]){
		printf("dBA[%d]=%d, B[%d]=%d. Program Stops!\n",i,dBA[i],i,B[i]);
		exit(1);
	}
	for(j=0;j<Ntp;j++)
		for(k=0;k<Modulation;k++){
			if(dBM[i][j][k]!=A[i][j][k]){
				printf("dBM[%d][%d][%d]=%d, A[%d][%d][%d]=%d. Program Stops!\n",i,j,k,dBM[i][j][k],i,j,k,A[i][j][k]);
				exit(1);
			}
	}	
}
double SE[Tbs]={0};
for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++)
		for(k=0;k<Modulation;k++)
			if(dBM[i][j][k]==1){
				CumulativeP[i]+=P[i][j][k];
				TotalP+=P[i][j][k];
				CumulativeDS[i]+=DataSub[j][k];
					CumulativeBR[i]+=BR[j];
				TotalDS+=DataSub[j][k];
					TotalBR+=BR[j];
				CumSerMS[i]++;
				TotSerMS++;
					NOV+=(BR[j]/DataSub[j][k]);
					SE[i]+=(BR[j]/DataSub[j][k]);
				if(dBA[i]!=1){
					printf("BS[%d] should be in active mode! Program stops.",i);
					exit(1);
				}	
			}
for(i=0;i<Tbs;i++){
	if(i<Tbs-mNbs)
		if(CumulativeP[i]>0){
			CumulativeP[i]+=ABP;
			TotalP+=ABP;
		}
		else{
			CumulativeP[i]+=SBP;
			TotalP+=SBP;
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(dBM[i][j][k]!=0){
						printf("BS%d should work in sleep mode! Program stops.",i);
						exit(1);
					}		
		}
	else
		if(CumulativeP[i]>0){
			CumulativeP[i]+=mABP;
			TotalP+=mABP;
		}
		else{
			CumulativeP[i]+=mSBP;
			TotalP+=mSBP;
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(dBM[i][j][k]!=0){
						printf("BS%d should work in sleep mode! Program stops.",i);
						exit(1);
					}		
		}
}

ConnectionPrinter(CplexOut);
//InfoPrinter(CplexOut);

for(i=0;i<Tbs;i++){
	if(i==0)
		fprintf(CplexOut,"Macro BS Statistics\n");
	else if(i==Tbs-mNbs)
		fprintf(CplexOut,"Micro BS Statistics\n");
       fprintf(CplexOut,"Power from BS %d: %lf\n",i,CumulativeP[i]);
       fprintf(CplexOut,"DS from BS %d: %d\n",i,CumulativeDS[i]);
       //if (i<Tbs-mNbs)
		fprintf(CplexOut, "Subframe fraction: %lf\n", x[i]/t);
	fprintf(CplexOut,"Bandwidth from BS %d: %lf\n",i,CumulativeBR[i]);
       fprintf(CplexOut,"MSs served by BS %d: %d\n",i,CumSerMS[i]);
	fprintf(CplexOut,"Spectrum efficiency: %lf\n",SE[i]);
	fprintf(CplexOut,"Spectrum-energy efficiency: %lf\n\n",SE[i]/CumulativeP[i]);
	
	//else
		//fprintf(CplexOut, "Subframe fraction: \lf\n", x[i]/t);
       }
  double TranPower = 0.0;
  for(i=0;i<Tbs;i++)
    for(j=0;j<Ntp;j++)
      for(k=0;k<Modulation;k++)
	if(dBM[i][j][k]==1)
	  TranPower += P[i][j][k];
  fprintf(CplexOut, "Transmit Power: %g\n", TranPower);
       fprintf(CplexOut,"Total Power %lf\n",TotalP);
       fprintf(CplexOut,"Total DSs: %d\n",TotalDS);
       fprintf(CplexOut,"Total served MSs: %d\n",TotSerMS);
	fprintf(CplexOut,"Total Bandwidth: %lf\n",TotalBR);
	fprintf(CplexOut,"Objective Value: %lf\n",TotalBR/TotalDS/TotalP);
	fprintf(CplexOut,"New Objective Value: %lf\n",NOV/TotalP);
	fprintf(CplexOut,"Transmit Objective Value: %lf\n",NOV/TranPower);
fclose(CplexOut);
//system("pause");
}

