#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"header/CoverageComputation.h"
#include<math.h>

#define Tbs 23
#define Ntp 4600
#define BP 0.08
#define SBP 8.00	//Power for BS in sleep mode
#define ABP 68.73	//Basic power for BS in active mode
#define Modulation 7
#define mNbs 8
#define mSBP 7.3	//Power for micro BSs in sleep mode
#define mABP 12	//Basic power for micro BSs in active mode
int Xb[Tbs]={0};
int Yb[Tbs]={0};
int Xtp[Ntp]={0};
int Ytp[Ntp]={0};
int dBM[Tbs][Ntp][Modulation];
int *ptrdBM;
int dBA[Tbs]={0};
double P[Tbs][Ntp][Modulation];
double BR[Ntp]={0};
int DataSub[Ntp][Modulation];
double CumulativeP[Tbs]={0};
double CumulativeBR[Tbs]={0};
int CumulativeDS[Tbs]={0};
int CumSerMS[Tbs]={0};
double TotalP=0.0;
double TotalBR=0.0;
double NOV=0.0;
int TotalDS=0;
int TotSerMS=0;
void ConnectionPrinter(FILE* CplexOut, int* ptrdBM);
/*
InfoPrinter(FILE* CplexOut){
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
}
*/
void ConnectionPrinter(FILE* CplexOut, int* ptrdBM){
               int i,j,k;
	double *coverage;
	
	coverage=CoverageComputation(ptrdBM,Tbs,Ntp,Xb,Yb,Xtp,Ytp,Modulation);

fprintf(CplexOut,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);

fprintf(CplexOut,"# plot \"PowerSavingCPLEX_Output.txt\" index 0:0 using 2:3:1 notitle with labels, \"PowerSavingCPLEX_Output.txt\" index 1:1 using 5:6:(($10)*100) title \"Energy demand of associated user\" with circles fill solid lt rgb \"red\", \"PowerSavingCPLEX_Output.txt]\" index 2:2 using 5:6:(($10)*100) title \"Energy demand of blocked user\" with circles fill solid lt rgb \"orange\"");//, \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
/*for(i=0;i<Tbs;i++)
fprintf(CplexOut,", \"PowerSavingCPLEX_Output.txt\" index %d:%d using 5:6:10 notitle with circles",i,i);*/
fprintf(CplexOut,"\n\n");

for(i=0;i<Tbs;i++)
	fprintf(CplexOut," \"BS[%d]\" %d %d %g\n",i,Xb[i],Yb[i],*(coverage+i));
fprintf(CplexOut,"\n\n");

	for(i=0;i<Tbs;i++){
		fprintf(CplexOut,"#BS[%d] %d %d BS[%d] %d %d\n",i,Xb[i],Yb[i],i,Xb[i],Yb[i]);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(dBM[i][j][k]==1)
					//fprintf(CplexOut,"BS[%d] %d %d BS[%d] %d %d\n",i,Xb[i],Yb[i],i,Xb[i],Yb[i]);
					fprintf(CplexOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW with modulation %d\n",i,Xb[i],Yb[i],j,Xtp[j],Ytp[j],DataSub[j][k],P[i][j][k],k);
				
			
	
		fprintf(CplexOut,"\n\n");
	}
}
                                  
main(){
       
FILE *infile, *outfile; 
int i,j,k;
char id[100];
char tempC;
int a, b, c;
int BS,MS,value,Mod;
size_t loc;
char string[12]="<variables>" ;
char string1[13]="</variables>" ;
char string2[10]="<variable" ;
	ptrdBM=&dBM[0][0][0];
if ((infile=fopen("PSCRD.txt", "r")) == NULL)
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
				}
			}
	while (fscanf(infile, "%s", &id)!=EOF){
		if (strcmp(id,string)==0){
			fscanf(infile, "%s", &id);
			while(strcmp(id,string1)!=0){
				//printf("%s\n",id);
				while(strcmp(id,string2)==0){//
					fscanf(infile, " name=\"");
					fscanf(infile, "%c", &tempC);
					//printf("%c\n",tempC);
					if(tempC=='X'){
						fscanf(infile, "%s",&id);
						//printf("%s\n",id);
						fscanf(infile, "%s",&id);
						//printf("%s\n",id);
						fscanf(infile, "%s",&id);
						//printf("%s\n",id);
					}
					else if(tempC=='d'){
						fscanf(infile, "%c%c", &tempC,&tempC);
						//printf("%c\n",tempC);
						if(tempC=='M'){
							fscanf(infile, "%d_%d_%d\"", &BS, &MS, &Mod);
							fscanf(infile, "%s", &id);
							fscanf(infile, " value=\"%d\"/>", &value);
							//printf("%d %d %d %d\n",BS,MS,Mod,value);
							dBM[BS][MS][Mod]=value;
							
						}
						else if(tempC=='A'){
							fscanf(infile, "%d\"", &BS);
							fscanf(infile, "%s", &id);
							fscanf(infile, " value=\"%d\"/>", &value);
							dBA[BS]=value;
						}
						else{
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
							fscanf(infile, "%s",&id);
							//printf("%s\n",id);
						}
					}
					fscanf(infile, "%s", &id);//printf("%s\n",id); 
					/*fscanf(infile, " name=\"dBM%d_%d\"", &BS, &MS);
					fscanf(infile, "%s", &id);
					fscanf(infile, " value=\"%d\"/>", &value);
					dBM[BS-1][MS-1]=value;
					fscanf(infile, "%s", &id);
					//printf("Inner Scanning\n");
					if(BS==Tbs&&MS==Ntp)
						strcpy(id,string1);*/
				}      //
				                                                    
			}    
			                                                       
		}//fscanf(infile, "%s", &id);printf("%s\n",id); 
	
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
char id1[100];

if ((CplexIn=fopen("outfile1.txt", "r")) == NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");
if ((CplexOut=fopen("PowerSavingCPLEX_Output.txt", "w")) == NULL)
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
      fscanf(CplexIn," with modulation %d", &Mod);
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
				if(dBA[i]!=1){
					printf("BS%d should operate in active mode! Program stops.",i);
					exit(1);
				}
			}
for(i=0;i<Tbs;i++){
	if(CumulativeP[i]>0&&i<Tbs-mNbs){
		CumulativeP[i]+=ABP;
		TotalP+=ABP;
	}
	else if(CumulativeP[i]>0&&i>=Tbs-mNbs){
		CumulativeP[i]+=mABP;
		TotalP+=mABP;
	}
	else if(CumulativeP[i]==0&&i<Tbs-mNbs){
		CumulativeP[i]+=SBP;
		TotalP+=SBP;
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(dBM[i][j][k]==1){
					printf("BS%d should run in sleep mode! Program stops.",i);
					exit(1);
				}
	}
	else if(CumulativeP[i]==0&&i>=Tbs-mNbs){
		CumulativeP[i]+=mSBP;
		TotalP+=mSBP;
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(dBM[i][j][k]==1){
					printf("BS%d should run in sleep mode! Program stops.",i);
					exit(1);
				}
	}
}

ConnectionPrinter(CplexOut,ptrdBM);
//InfoPrinter(CplexOut);

	for(i=0;i<Tbs;i++){
		if(i==0)
			fprintf(CplexOut,"Macro BS Statistics\n");
		if(i==Tbs-mNbs)
			fprintf(CplexOut,"Micro BS Statistics\n");
	fprintf(CplexOut,"Power from BS %d: %lf\n",i,CumulativeP[i]);
	fprintf(CplexOut,"DS from BS %d: %d\n",i,CumulativeDS[i]);
		fprintf(CplexOut,"Bandwidth from BS %d: %lf\n",i,CumulativeBR[i]);
	fprintf(CplexOut,"MSs served by BS %d: %d\n\n",i,CumSerMS[i]); 
	}
       fprintf(CplexOut,"Total Power %lf\n",TotalP);
       fprintf(CplexOut,"Total DSs: %d\n",TotalDS);
       fprintf(CplexOut,"Total served MSs: %d\n",TotSerMS);
	fprintf(CplexOut,"Total Bandwidth: %lf\n",TotalBR);
	fprintf(CplexOut,"Objective Value: %lf\n",TotalBR/TotalDS/TotalP);
	fprintf(CplexOut,"New Objective Value: %lf\n",NOV/TotalP);
fclose(CplexOut);
//system("pause");
}

