//#ifndef MY_HEADER_H
//#define MY_HEADER_H
//extern int Ntp;
//extern int Modulation;
//extern int mNbs;
//extern int DSt;
extern int* ptrAdjacentBSs;
extern int* ptrDS;
extern int* ptrminM;

//extern double mABP;
//extern double mSBP;
//extern double ABP;
//extern double SBP;
//extern double MP;
//extern double mP;
//extern double BP;
extern double* ptrBR;
extern double* ptrP;

localOptimal(int Tbs, int nUser, int Ntp, int Modulation, int mNbs, int DSt, double mSBP, double mABP, double ABP, double SBP, double MP, double mP, double BP){
  int i, j, l;
  int* connection = malloc( Ntp * sizeof(int));
  int* operation = malloc( Tbs * sizeof(int));
  int* association = malloc( Tbs * Ntp * Modulation * sizeof(int));
  //int totalCalls = 0;
  printf("Start the local optimal huristic...\n");
 // memset(connection, 0, Ntp);
 // memset(operation, 0, Tbs);
 // memset(association, 0, Tbs * Ntp * Modulation);
  clock_t start_tick_all, end_tick_all;
  
  for (i = 0; i < Tbs; i++){
    *(operation + i) = 0;
    if (*(operation + i) != 0){
      printf("wtf(%d)\n", *(operation + i));
      exit(1);
    }
    for (j = 0; j < Ntp; j++){
      *(connection + j) = 0;
      if (*(connection + j) != 0){
	printf("wtf(%d, %d)\n", j, *(connection + j));
	exit(1);
      }
      for (l = 0; l < Modulation; l++){
	*(association + i * Ntp * Modulation + j * Modulation + l) = 0;
	if (*(association + i * Ntp * Modulation + j * Modulation + l) != 0){
	  printf("wtf\n");
	  exit(1);
	}
      }
    }
  }
  start_tick_all = clock();
  //omp_set_nested(1);
  #pragma omp parallel for default(none) private(j, l) schedule(dynamic) shared(Tbs, mNbs, ptrAdjacentBSs, Ntp, connection, mP,DSt, Modulation, mABP, operation, association, ptrDS)
  for (i = Tbs - mNbs; i < Tbs; i++){
    int UENum = 0;
    clock_t start_tick, end_tick;
    double elapsed;
    //int ID = omp_get_thread_num();   
    start_tick=clock();
    for (j = 0; j < Ntp; j++)
      if (*(ptrAdjacentBSs + i * Ntp + j) == 1 && *(connection + j) == 0)
	UENum++;
    if (UENum > 0){
      double currentOpt = 0.0;
      int* currentUE = NULL;
      int* currentMod = NULL;
      int currentUENum = 0;
      printf("\neNodeB%d has %d UEs in coverage(%d).\n", i, UENum, omp_get_thread_num());
      //#pragma omp parallel for default(none) schedule(static) shared(UENum, i, mP, DSt, Modulation, mABP, Ntp, connection, currentUENum, currentMod, currentUE, currentOpt)
      //for (l = 1; l <= UENum; l++){
      l = UENum;
	int callingTimes = 0;
	printf("Working on the %d UEs of eNodeB %d (%d).\n", l, i, omp_get_thread_num());
	SEE_UE(l, &currentOpt, &currentUE, &currentMod, &currentUENum, i, &callingTimes, connection, Ntp, mABP, Modulation, DSt, mP);
	time_t t = time(NULL);
	struct tm *tm = localtime(&t);
	printf("\nWorked on %d of %d UEs of eNodeB %d for (%d)asctime = %s", l, UENum, i, omp_get_thread_num(), asctime(tm));/*
	int numinator = UENum;
	int numinatorVar = UENum;
	int m = l;
	int denuminator = 1;
	int denuminatorVar = 1;
	while (m > 1){
	  numinator = numinator * (--numinatorVar);
	  denuminator = denuminator * (++denuminatorVar);
	  m--;
	}
	if (callingTimes != (numinator / denuminator) * (int) pow(Modulation, l)){
	  printf("calling times wrong (%d, %d, %d, %d, %d, %d)\n", l, UENum, callingTimes, (numinator) , (denuminator), (int) pow(Modulation, l));
	  exit(1);
	}*/
      //}
      if (currentUENum == 0 || currentUENum > UENum || currentUE == NULL || (!currentMod) || currentOpt == 0.0){
	printf("cannot do this trick\n");
	exit(1);
      }
      int RBTotal = 0;
      for (l = 0; l < currentUENum; l++){
	*(operation + i) = 1;
	*(association + i * Ntp * Modulation + *(currentUE + l) * Modulation + *(currentMod + l)) = 1;
	*(connection + *(currentUE + l)) = 1;
	//printf("%d %d %d\n", i, *(currentUE + l), *(currentMod + l));
	RBTotal += *(ptrDS + *(currentUE + l) * Modulation + *(currentMod + l));
      }
     // printf("%d %d consumes %d RBs.\n", i, *(operation + i), RBTotal);
      printf("%d serves %d UEs(%d).\n", i, currentUENum, omp_get_thread_num());
      free(currentUE);
      free(currentMod);
    }
    //printf("%d calls for eNodeB%d\n", callingTimes, i);
  //  totalCalls += callingTimes;
    end_tick=clock();
    elapsed=(double)(end_tick-start_tick)/CLOCKS_PER_SEC;
    printf("Running Time for Dealing with femto%d: %g secs(by thread%d)\n", i, elapsed, omp_get_thread_num());
  }
  end_tick_all = clock();
  printf("Running Time for all the femtos: %g mins\n\n", (double)(end_tick_all-start_tick_all)/CLOCKS_PER_SEC/ 60.0);
  locOptModel(Tbs, nUser, operation, association, Ntp, Modulation, mNbs, DSt, MP, mP, BP, ABP, SBP, mABP, mSBP);
  free(connection);
  free(operation);
  free(association);
}

SEE_UE(int num, double* currentOpt, int** currentUE, int** currentMod, int* currentUENum, int eNodeB, int* callingTimes, int* connection, int Ntp, double mABP, int Modulation, int DSt, double mP){
  int i;
  int* UEArr = NULL;
  int end = 0;
  do{
    end = update(&UEArr, num, eNodeB, connection, Ntp);
    //for (i = 0; i < num; i++)
      //printf("%d\n", *(UEArr + i));
    if (UEArr == NULL){
      
      printf("why(%d)?\n", *(UEArr + 0));
      exit(1);
    }
    if (end == 0){
      SEE_Mod(UEArr, currentOpt, currentUE, currentMod, currentUENum, eNodeB, num, callingTimes, mABP, Modulation, Ntp, DSt, mP);
      //printf("%d\t", *callingTimes);
    }
  }while (end == 0);
  free(UEArr);
}

int update(int** UEArr, int size, int eNodeB, int* connection, int Ntp){
  int j;
  int end = 0;
  if (*UEArr == NULL){
    //printf("constructing.....\n");
    *UEArr = malloc(sizeof(int) * size);
    //printf("%x\n", *UEArr);
    int ind =0;
    for (j = 0; j < Ntp; j++){
      if (*(ptrAdjacentBSs + eNodeB * Ntp + j) == 1 && *(connection + j) == 0){
	if (ind < size){
	  *(*UEArr + ind) = j;
	  ind++;
	  //printf("%d\t", **(UEArr + (ind - 1)));
	  //printf("test: %d\t", *(*UEArr + (ind - 1)));
	}
	
      }
    }
  }
  else
    end = findNextUE(*UEArr, size, size - 1, eNodeB, connection, Ntp);
  
  return end;
}

int findNextUE(int* UEArr, int size, int ind, int eNodeB, int* connection, int Ntp){
  int j;
  int nextUE = -1;
  int end = 0;
  if (ind < 0){
    printf("index error\n");
    exit(1);
  }
  for (j = *(UEArr + ind) + 1; j < Ntp; j++){
    //printf("%d\t", j);
    if (j < Ntp)
      if (*(ptrAdjacentBSs + eNodeB * Ntp + j) == 1 && *(connection + j) == 0)
	if (ind < size - 1){
	  if (j < *(UEArr + (ind + 1)))
	    nextUE = j;
	  break;
	}
	else{
	  nextUE = j;
	  break;
	}
  }
  if (nextUE != -1){
    *(UEArr + ind) = nextUE;
    if (ind < size - 1)
      updateOrd(UEArr, size, ind, eNodeB, connection, Ntp);
  }
  else{
    
    if (ind == 0)
      end = 1;
    else{
      //printf("%d\n", ind);
      end = findNextUE(UEArr, size, ind - 1, eNodeB, connection, Ntp);
    }
  }
  return end;
}

updateOrd(int* UEArr, int size, int ind, int eNodeB, int* connection, int Ntp){
  int pos = ind + 1;
  int j;
  for (j = *(UEArr + ind) + 1; j < Ntp; j++){
    if (*(ptrAdjacentBSs + eNodeB * Ntp + j) == 1 && *(connection + j) == 0){
      *(UEArr + pos) = j;
      pos++;
    }
    if (pos == size)
      break;
  }  
}

double calculateSEE(int* ModArr, int* UEArr, int size, int eNodeB, double mABP, int Modulation, int Ntp, int DSt, double mP){
  int l;
  double NOV = 0.0;
  double power = mABP;
  int totalRB = 0;
  //double totalPower = 0.0;
  for (l = 0; l < size; l++){
    totalRB += *(ptrDS + *(UEArr + l) * Modulation + *(ModArr + l));
    //totalPower += *(ptrP + eNodeB * Ntp * Modulation + *(UEArr + l) * Modulation + *(ModArr + l));
    NOV += (*(ptrBR + *(UEArr + l))) / (*(ptrDS + *(UEArr + l) * Modulation + *(ModArr + l)));
    power += *(ptrP + eNodeB * Ntp * Modulation + *(UEArr + l) * Modulation + *(ModArr + l));
  }
  if (totalRB > DSt || (power - mABP) > mP)
    return 0.0;
  else
    return NOV / power;
}

SEE_Mod(int* UEArr, double* currentOpt, int** currentUE, int** currentMod, int* currentUENum, int eNodeB, int size, int* callingTimes, double mABP, int Modulation, int Ntp, int DSt, double mP){
  double SEE = 0.0;
  int* ModArr = NULL;
  int end = 0;
  do{
    end = nextMod(&ModArr, size, callingTimes, Modulation);
   // printf("\n%d\t", *callingTimes);
    if (ModArr == NULL){
      printf("why?\n");
      exit(1);
    }
    if (end == 0){
      SEE = calculateSEE(ModArr, UEArr, size, eNodeB, mABP, Modulation, Ntp, DSt, mP);
      compareSEE(SEE, currentOpt, currentUE, currentMod, ModArr, UEArr, size, currentUENum);
    }
  }while(end == 0);
  free(ModArr);
}

int nextMod(int** ModArr, int size, int* callingTimes, int Modulation){
  int l;
  int end = 0;
  if (*ModArr == NULL){
    *ModArr = malloc(sizeof(int) * size);
    for (l = 0; l < size; l++)
      *(*ModArr + l) = 0;
    *callingTimes = *callingTimes + 1;
  }
  else
    end = updateMod(*ModArr, size, size - 1, callingTimes, Modulation);
  return end;
}

int updateMod(int* ModArr, int size, int ind, int* callingTimes, int Modulation){
  int end = 0;
  if (*(ModArr + ind) == Modulation - 1)
    if (ind != 0)
      end = updateMod(ModArr, size, ind - 1, callingTimes, Modulation);
    else
      end = 1;
  else{
    *(ModArr + ind) = *(ModArr + ind) + 1;
    shuffleMod(ModArr, ind, size);
    *callingTimes = *callingTimes + 1;
  }
  return end;
}

shuffleMod(int* ModArr, int ind, int size){
  int l;
  if (ind < size - 1)
    for (l = ind + 1; l < size; l++)
      *(ModArr + l) = 0;
}

double my_round(double x, unsigned int digits){
  if (digits > 0){
    return my_round(x*10.0, digits-1)/10.0;
  }
  else{
    return floor(x);
  }
}

compareSEE(double SEE, double* currentOpt, int** currentUE, int** currentMod, int* ModArr, int* UEArr, int size, int* currentUENum){
  int l;
  if (SEE > *currentOpt){
    #pragma omp critical
    if (SEE > *currentOpt){/*
    printf("previous SEE %g, new SEE %g\n", *currentOpt, SEE);
    for (l = 0; l < *currentUENum; l++){
      printf("(%d, %d)\t", *(*currentMod + l), *(*currentUE + l));
    }
    printf("\n");*/
    
      *currentOpt = SEE;
      *currentUENum = size;
      *currentMod = realloc(*currentMod, sizeof(int) * size);
      *currentUE = realloc(*currentUE, sizeof(int) * size);
      for (l = 0; l < size; l++){
	*(*currentMod + l) = *(ModArr + l);
	*(*currentUE + l) = *(UEArr + l);
      }/*
    for (l = 0; l < size; l++){
      printf("(%d, %d)\t", *(*currentMod + l), *(*currentUE + l));
    }
    printf("\n");*/
    }
  }/*
  else{
    printf("previous SEE %g, calculated SEE %g\n", *currentOpt, SEE);
    for (l = 0; l < *currentUENum; l++){
      printf("(%d, %d)\t", *(*currentMod + l), *(*currentUE + l));
    }
    printf("\n");
    for (l = 0; l < size; l++){
      printf("(%d, %d)\t", *(ModArr + l), *(UEArr + l));
    }
    printf("\n");
  }*/
}

locOptModel(int Tbs, int nUser, int* operation, int* association, int Ntp, int Modulation, int mNbs, int DSt, double MP, double mP, double BP, double ABP, double SBP, double mABP, double mSBP){
	FILE *PrintToFile;
	if((PrintToFile=fopen("LocOptModel.txt","w"))==NULL)
		printf("\nerror!Fail to open LocOptModel.txt!");
	else
		printf("\nOpen LocOptModel.txt successfully!\n");
	fprintf(PrintToFile,"This is the input to CPLEX for Optimization Model with local optimal heuristic.\n");
	
	cobjective(Tbs,PrintToFile,Ntp,Modulation,ptrBR,ptrDS);
  int i, j, k, l;
  for (i = Tbs- mNbs; i < Tbs; i++){
    //printf("%d, %d\n", i, *(operation + i));
    int associationNumber = 0;
    for (j = 0; j< Ntp; j++){
      int associationConstraint = 0;
      for (k = 0; k < Modulation; k++){
	fprintf(PrintToFile, "a%d_%d_%d = %d\n", i, j, k, *(association + i * Ntp * Modulation + j * Modulation + k));
	associationConstraint += *(association + i * Ntp * Modulation + j * Modulation + k);
	associationNumber += *(association + i * Ntp * Modulation + j * Modulation + k);
      }
      if (associationConstraint > 1){
	printf("association constraint error\n");
	exit(1);
      }
    }
    if ((*(operation + i) == 0 && associationNumber != 0) || (*(operation + i) != 0 && associationNumber == 0)){
      printf("association number error\n");
      exit(1);
    }
  }/*
  for (i = Tbs - mNbs; i < Tbs; i++){
   // fprintf(PrintToFile, "b%d = %d\n", i, *(operation + i));
    printf("%d, %d\n", i, *(operation + i));
  }*/
  for (i = Tbs - mNbs; i < Tbs; i++){
    fprintf(PrintToFile, "b%d = %d\n", i, *(operation + i));
    //printf("%d, %d\n", i, *(operation + i));
  }
  for (i = 0; i < Tbs - mNbs; i++){
    int highestRBNum = 0;
    for (l = Tbs - mNbs; l < Tbs; l++)
      if (*(ptrminM + i * Tbs + l) == 1){
	int RBNum = 0;
	for (j = 0; j < Ntp; j++)
	  for (k = 0; k < Modulation; k++)
	    if (*(association + l * Ntp * Modulation + j * Modulation + k) == 1)
	      RBNum += *(ptrDS + j * Modulation + k);
	if (RBNum > highestRBNum)
	  highestRBNum = RBNum;
      }
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, my_round(((double)(DSt - highestRBNum) / (double)DSt), 5));
    //fprintf(PrintToFile, "x%d - %lf t = 0\n", i, ((double)(DSt - highestRBNum) / (double)DSt));
    //printf("%lf\t", ((double)(DSt - highestRBNum) / (double)DSt));
  }
	
	ABSconstraint1(Tbs,ptrDS,PrintToFile,ptrminM,mNbs,Modulation,ptrDS,DSt,Ntp);
	
	cconstraint2(Tbs,PrintToFile, mNbs, Ntp, Modulation, ptrP, MP, mP);
	cconstraint3(Tbs,PrintToFile,nUser, Modulation);
	cconstraint4(Tbs,PrintToFile, Ntp, Modulation, BP);
	cconstraint5(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint10(Tbs,PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
	cconstraint6(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint7(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint8(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint9(Tbs,PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, ptrP, mABP, mSBP);
	cconstraint11(Tbs,PrintToFile,nUser, Ntp, Modulation);

	ABSbounds(Tbs,PrintToFile, Ntp, Modulation, mNbs);
	
	cspecifyTypes(Tbs,PrintToFile, Ntp, Modulation);
	fclose(PrintToFile);
}  


//#endif 
