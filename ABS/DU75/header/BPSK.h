
void BPSK(int Tbs, int* DS, double* P, double* BR, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int* Connection, int Ntp, double BP, int Modulation, int MP, int DSt, double mABP, double mSBP, int mNbs, double mp, int* minM, int Rbs, int rbs, int space, double ABP, double SBP){
	int i,j,k=0;
	double TotalBandwidth=0.0, TotalPower=0.0, NOV=0.0;
	int TotalDS=0;
	double *Power=malloc(Tbs*sizeof(double));
	int *DataSub=malloc(Tbs*sizeof(int));
	for(i=0;i<Tbs;i++){
		*(Power+i)=0;
		*(DataSub+i)=0;
	}
	FILE *BPSK;
	if((BPSK=fopen("BPSK.txt","w"))==NULL){
		printf("\n\nerror!Fail to open file!");
		exit(1);
	}
	else
		printf("\n\nOpen BPSK successfully!\n");
	fprintf(BPSK,"\n#This file expresses fixed cell sizes with BPSK.\n");
	fprintf(BPSK,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	for(i=0;i<Tbs-mNbs;i++)
		fprintf(BPSK," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],Rbs/space);
	for(i=Tbs-mNbs;i<Tbs;i++)
		fprintf(BPSK," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],rbs/space);
	fprintf(BPSK,"\n\n");
	for(i=0;i<Tbs;i++){
		fprintf(BPSK,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
		for(j=0;j<Ntp;j++)
			if(*(Connection+j)<0){
				printf("MS has no connections! Program Stops!!\n");
				exit(1);
			}
			else if(*(Connection+j)==i){
				Power[i]+=*(P+i*Ntp*Modulation+j*Modulation+k);
					if(i<Tbs-mNbs&&Power[i]>MP){
						printf("BS%d has too much power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					else if(i>=Tbs-mNbs&&Power[i]>mp){
						printf("BS%d has too much power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
				DataSub[i]+=*(DS+j*Modulation+k);
					if(DataSub[i]>DSt){
						printf("BS%d has too many subcarrier demands! Program Stops!\n",i);
						exit(1);
					}
				TotalBandwidth+=BR[j];
				TotalDS+=*(DS+j*Modulation+k);
				TotalPower+=*(P+i*Ntp*Modulation+j*Modulation+k);
				
				NOV+=(*(BR+j)/(*(DS+j*Modulation+k)));
				fprintf(BPSK,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
				fprintf(BPSK,"BS[%d] %d %d MS[%d] %d %d DS %d power %lf  mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),BR[j],k);
			}
//		if(Power[i]>0){
			if(i<Tbs-mNbs){
				Power[i]+=ABP;
				TotalPower+=ABP;
			}
			else{
				Power[i]+=mABP;
				TotalPower+=mABP;
			}
//		}
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
	int temp;
	for(i=0;i<Tbs-mNbs;i++){
		temp=DataSub[i];
		for(j=Tbs-mNbs;j<Tbs;j++)
			if(*(minM+i*Tbs+j)==1)
				temp+=DataSub[j];
		if(temp>DSt){
			printf("The group of BS%d has too many subcarrier demands. Program Stops!\n",i);
			exit(1);
		}				
	}
	fprintf(BPSK,"\n\n");
	fprintf(BPSK,"Total Power: %lf\n",TotalPower);
	fprintf(BPSK,"Total DSs: %d\n",TotalDS);
	fprintf(BPSK,"Total BR: %lf\n",TotalBandwidth);
	fprintf(BPSK,"Spectrum-energy Efficiency: %lf\n", NOV/TotalPower);
	fprintf(BPSK,"\n\n");
	fprintf(BPSK," plot \"BPSK.txt\" index 0:0 using 2:3:1 notitle with labels, \"BPSK.txt\" index 0:0 using 2:3:4 title \"users\" with points");
	for(i=1;i<=Tbs;i++)
		fprintf(BPSK,", \"BPSK.txt\" index %d:%d using 5:6 title \" BS[%d]\" with lines",i,i,i-1);
	fclose(BPSK);
	free(Power);
	free(DataSub);
}
