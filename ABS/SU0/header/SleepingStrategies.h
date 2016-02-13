extern double* ptrBR;
extern double* ptrP;
extern int* ptrAdjacentBSs;
extern int* ptrDS;
extern int* ptrminM;


int bernoulli(float p){
  if(p < 0 || p > 1) return -1;
  float x = (float)rand()/(float)(RAND_MAX/1);
  if(p < x) return 0;
  return 1;
}

binary_sleeping(int Tbs, int nUser, int Ntp, int Modulation, int mNbs, int DSt, double mSBP, double mABP, double ABP, double SBP, double MP, double mP, double BP, double ABS){
  int i, j;
  int associated_UE_Num = 0;
  int* associated_UE = malloc(Ntp * sizeof(int));
  int* active = malloc(Tbs * sizeof(int));
  float binary_probability = 0.5;
  for(i = 0; i < Ntp; i++)
    *(associated_UE + i) = 0;
  for(i = 0; i < Tbs; i++)
    *(active + i) = 0;
  while(associated_UE_Num < Ntp * (1 - BP)){
    i = rand() % Tbs;
    if(*(active + i) == 0){
      *(active + i) = bernoulli(binary_probability);
      if(*(active + i) == 1)
	for(j = 0; j < Ntp; j++)
	  if(*(associated_UE + j) == 0 && *(ptrAdjacentBSs + i * Ntp + j) == 1){
	    *(associated_UE + j) = 1;
	    associated_UE_Num++;
	  }
      else if(*(active + i) == -1){
	printf("operation mode error.\n");
	exit(1);
      }
    }
  }
  free(associated_UE);
  char name[80];
  char fileName[80];
  char str[10];
  sprintf(str, "%g", ABS);
  strcpy(name, "Binary Sleeping(");
  strcat(name, str);
  strcat(name, ")");
  strcpy(fileName, "BinarySleeping(");
  strcat(fileName, str);
  strcat(fileName, ").txt");
  FILE *PrintToFile;
  if((PrintToFile = fopen(fileName, "w")) == NULL)
    printf("\nerror! Fail to open %s!", fileName);
  else
    printf("\nOpen %s successfully!\n", fileName);
  fprintf(PrintToFile, "This is the input to CPLEX for %s.\n", name);
  cobjective(Tbs, PrintToFile, Ntp, Modulation, ptrBR, ptrDS);
  for(i = 0; i < Tbs - mNbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, 1 - ABS);
  for(i = Tbs - mNbs; i < Tbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, ABS);
  ABSconstraint1(Tbs, ptrDS, PrintToFile, ptrminM, mNbs, Modulation, ptrDS, DSt, Ntp);
  cconstraint2(Tbs, PrintToFile, mNbs, Ntp, Modulation, ptrP, MP, mP);
  cconstraint3(Tbs, PrintToFile, nUser, Modulation);
  cconstraint4(Tbs, PrintToFile, Ntp, Modulation, BP);
  //cconstraint5(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint6(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint7(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint8(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint9(Tbs, PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, ptrP, mABP, mSBP);
  cconstraint10(Tbs, PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
  cconstraint11(Tbs, PrintToFile, nUser, Ntp, Modulation);
  for(i = 0; i < Tbs; i++)
    fprintf(PrintToFile, "b%d = %d\n", i, *(active +i));
  ABSbounds(Tbs, PrintToFile, Ntp, Modulation, mNbs);
  cspecifyTypes(Tbs, PrintToFile, Ntp, Modulation);
  fclose(PrintToFile);
  free(active);
} 


uniform_sleeping(int Tbs, int nUser, int Ntp, int Modulation, int mNbs, int DSt, double mSBP, double mABP, double ABP, double SBP, double MP, double mP, double BP, double ABS){
  int i, j;
  int associated_UE_Num = 0;
  int* associated_UE = malloc(Ntp * sizeof(int));
  int* active = malloc(Tbs * sizeof(int));
  float* binary_probability = malloc(Tbs * sizeof(float));
  int* RB = malloc(Tbs * sizeof(int));
  for(i = 0; i < Ntp; i++)
    *(associated_UE + i) = 0;
  for(i = 0; i < Tbs; i++){
    *(active + i) = 0;
    *(RB + i) = 0;
  }
  for(i = 0; i < Tbs; i++)
    *(binary_probability + i) = ((float)(rand() % 10) + 1) / 10.0;
  while(associated_UE_Num < Ntp * (1 - BP)){
    i = rand() % Tbs;
    if(*(active + i) == 0){
      *(active + i) = bernoulli(*(binary_probability + i));
      if(*(active + i) == 1)
	for(j = 0; j < Ntp; j++)
	  if(*(associated_UE + j) == 0 && *(ptrAdjacentBSs + i * Ntp + j) == 1){
	    *(associated_UE + j) = 1;
	    *(RB + i) = *(ptrDS + j * Modulation + 0);
	    associated_UE_Num++;
	  }
      else if(*(active + i) == -1){
	printf("operation mode error.\n");
	exit(1);
      }
    }
  }
  for(i = 0; i < Tbs; i++)
    if(i < Tbs - mNbs && *(RB + i) > DSt * (1 - ABS)){
      printf("RB error\n");
      exit(1);
    }
    else if(i >= Tbs - mNbs && *(RB + i) > DSt *  ABS){
      printf("RB error\n");
      exit(1);
    }
  free(RB);
  free(associated_UE);
  free(binary_probability);
  char name[80];
  char fileName[80];
  char str[10];
  sprintf(str, "%g", ABS);
  strcpy(name, "Uniform Sleeping(");
  strcat(name, str);
  strcat(name, ")");
  strcpy(fileName, "UniformSleeping(");
  strcat(fileName, str);
  strcat(fileName, ").txt");
  FILE *PrintToFile;
  if((PrintToFile = fopen(fileName, "w")) == NULL)
    printf("\nerror! Fail to open %s!", fileName);
  else
    printf("\nOpen %s successfully!\n", fileName);
  fprintf(PrintToFile, "This is the input to CPLEX for %s.\n", name);
  cobjective(Tbs, PrintToFile, Ntp, Modulation, ptrBR, ptrDS);
  for(i = 0; i < Tbs - mNbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, 1 - ABS);
  for(i = Tbs - mNbs; i < Tbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, ABS);
  ABSconstraint1(Tbs, ptrDS, PrintToFile, ptrminM, mNbs, Modulation, ptrDS, DSt, Ntp);
  cconstraint2(Tbs, PrintToFile, mNbs, Ntp, Modulation, ptrP, MP, mP);
  cconstraint3(Tbs, PrintToFile, nUser, Modulation);
  cconstraint4(Tbs, PrintToFile, Ntp, Modulation, BP);
  //cconstraint5(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint6(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint7(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint8(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint9(Tbs, PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, ptrP, mABP, mSBP);
  cconstraint10(Tbs, PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
  cconstraint11(Tbs, PrintToFile, nUser, Ntp, Modulation);
  for(i = 0; i < Tbs; i++)
    fprintf(PrintToFile, "b%d = %d\n", i, *(active + i));
  ABSbounds(Tbs, PrintToFile, Ntp, Modulation, mNbs);
  cspecifyTypes(Tbs, PrintToFile, Ntp, Modulation);
  fclose(PrintToFile);
  free(active);
}

dynamic_activity_criterion_sleeping(int Tbs, int nUser, int Ntp, int Modulation, int mNbs, int DSt, double mSBP, double mABP, double ABP, double SBP, double MP, double mP, double BP, double ABS){
  int i, j;
  int associated_UE_Num = 0;
  int* associated_UE = malloc(Ntp * sizeof(int));
  int* active = malloc(Tbs * sizeof(int));
  float* binary_probability = malloc(Tbs * sizeof(float));
  int* load_demand = malloc(Tbs * sizeof(int));
  for(i = 0; i < Ntp; i++)
    *(associated_UE + i) = 0;
  for(i = 0; i < Tbs; i++){
    *(active + i) = 0;
    *(load_demand + i) = 0;
  }
  for(i = 0; i < Tbs; i++)
    for(j = 0; j < Ntp; j++)
      if(*(ptrAdjacentBSs + i * Ntp + j) == 1)
	*(load_demand + i) = *(load_demand + i) + *(ptrDS + j * Modulation + 0);
  for(i = 0; i < Tbs; i++)
    if (i < Tbs - mNbs){
      *(binary_probability + i) = (float)*(load_demand + i) / (float)(DSt * (1 - ABS));
      if(*(binary_probability + i) > 1 || *(binary_probability + i) < 0){
	printf("eNodeB%d load demand error (%d, %g, %d, %g, %g)\n", i, *(load_demand + i), DSt * (1 - ABS), DSt, ABS, (1 - ABS));
	exit(1);
      }
    }
    else{
      *(binary_probability + i) = (float)*(load_demand + i) / (float)(DSt * ABS);
      if(*(binary_probability + i) > 1 || *(binary_probability + i) < 0){
	printf("eNodeB%d load demand error (%d, %g, %d, %g)\n", i, *(load_demand + i), DSt * ABS, DSt, ABS);
	exit(1);
      }
    }
  while(associated_UE_Num < Ntp * (1 - BP)){
    i = rand() % Tbs;/*
    int k;
    for(k = 0; k < Tbs; k++)
      printf("%d ", *(active + k));;
    printf("\teNodeB%d is seleced(%d).\n", i, associated_UE_Num);*/
    if(*(active + i) == 0){
      *(active + i) = bernoulli(*(binary_probability + i));
      //printf("%d is still sleeping with probability %g.\n", i, *(binary_probability + i));
      if(*(active + i) == 1)
	for(j = 0; j < Ntp; j++)
	  if(*(associated_UE + j) == 0 && *(ptrAdjacentBSs + i * Ntp + j) == 1){
	    *(associated_UE + j) = 1;
	    associated_UE_Num++;
	  }
      else if(*(active + i) == -1){
	printf("operation mode error.\n");
	exit(1);
      }
    }
  }
  free(associated_UE);
  free(binary_probability);
  free(load_demand);
  char name[80];
  char fileName[80];
  char str[10];
  sprintf(str, "%g", ABS);
  strcpy(name, "Dynamic Activity Level Sleeping(");
  strcat(name, str);
  strcat(name, ")");
  strcpy(fileName, "DynamicSleeping(");
  strcat(fileName, str);
  strcat(fileName, ").txt");
  FILE *PrintToFile;
  if((PrintToFile = fopen(fileName, "w")) == NULL)
    printf("\nerror! Fail to open %s!", fileName);
  else
    printf("\nOpen %s successfully!\n", fileName);
  fprintf(PrintToFile, "This is the input to CPLEX for %s.\n", name);
  cobjective(Tbs, PrintToFile, Ntp, Modulation, ptrBR, ptrDS);
  for(i = 0; i < Tbs - mNbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, 1 - ABS);
  for(i = Tbs - mNbs; i < Tbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, ABS);
  ABSconstraint1(Tbs, ptrDS, PrintToFile, ptrminM, mNbs, Modulation, ptrDS, DSt, Ntp);
  cconstraint2(Tbs, PrintToFile, mNbs, Ntp, Modulation, ptrP, MP, mP);
  cconstraint3(Tbs, PrintToFile, nUser, Modulation);
  cconstraint4(Tbs, PrintToFile, Ntp, Modulation, BP);
  cconstraint5(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint6(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint7(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint8(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint9(Tbs, PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, ptrP, mABP, mSBP);
  cconstraint10(Tbs, PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
  cconstraint11(Tbs, PrintToFile, nUser, Ntp, Modulation);
  for(i = 0; i < Tbs; i++)
    fprintf(PrintToFile, "b%d = %d\n", i, *(active + i));
  ABSbounds(Tbs, PrintToFile, Ntp, Modulation, mNbs);
  cspecifyTypes(Tbs, PrintToFile, Ntp, Modulation);
  fclose(PrintToFile);
  free(active);
}

static_activity_criterion_sleeping(int Tbs, int nUser, int Ntp, int Modulation, int mNbs, int DSt, double mSBP, double mABP, double ABP, double SBP, double MP, double mP, double BP, double ABS){
  int i, j;
  int associated_UE_Num = 0;
  int* associated_UE = malloc(Ntp * sizeof(int));
  int* active = malloc(Tbs * sizeof(int));
  float* binary_probability = malloc(Tbs * sizeof(float));
  int* load_demand = malloc(Tbs * sizeof(int));
  float activity_criterion = 1.0 / 4.0;
  float low_activity_prob = 1.0 / 2.0;
  double* power = malloc(Tbs * sizeof(double));
  for(i = 0; i < Ntp; i++)
    *(associated_UE + i) = 0;
  for(i = 0; i < Tbs; i++){
    *(active + i) = 0;
    *(load_demand + i) = 0;
  }
  for(i = 0; i < Tbs; i++)
    for(j = 0; j < Ntp; j++)
      if(*(ptrAdjacentBSs + i * Ntp + j) == 1){
	*(load_demand + i) = *(load_demand + i) + *(ptrDS + j * Modulation + 0);
	*(power + i) = *(ptrP + i * Ntp * Modulation + j * Modulation + 0);
      }
  for(i = 0; i < Tbs; i++)
    if(i < Tbs - mNbs)
      if(*(power + i) > MP){
	printf("power error");
	exit(1);
      }
    else
      if(*(power + i) > mP){
	printf("power error");
	exit(1);
      }
  for(i = 0; i < Tbs; i++)
    if(i < Tbs - mNbs){
      if((float)*(load_demand + i) > (float)(DSt * (1 - ABS)) || (float)*(load_demand + i) / (float)(DSt * (1 - ABS)) < 0.0){
	printf("eNodeB%d load demand error(%d, %g).\n", i, *(load_demand + i), DSt * (1 - ABS));
	exit(1);
      }
      if(*(load_demand + i) == 0)
	*(binary_probability + i) = 0.0;
      else if((float)*(load_demand + i) / (float)(DSt * (1 - ABS)) < activity_criterion)
	*(binary_probability + i) = low_activity_prob;
      else
	*(binary_probability + i) = 1;
    }
    else{
      if((float)*(load_demand + i) > (float)(DSt * ABS) || (float)*(load_demand + i) / (float)(DSt * ABS) < 0.0){
	printf("eNodeB%d load demand error(%d, %g).\n", i, *(load_demand + i), DSt * ABS);
	exit(1);
      }
      if(*(load_demand + i) == 0)
	*(binary_probability + i) = 0.0;
      else if((float)*(load_demand + i) / (float)(DSt * ABS) < activity_criterion)
	*(binary_probability + i) = low_activity_prob;
      else
	*(binary_probability + i) = 1;
    }
  while((float)associated_UE_Num < (float)Ntp * (1 - BP)){
    i = rand() % Tbs;/*
    int k;
    for(k = 0; k < Ntp; k++)
      if(*(associated_UE + k) == 1)
	printf("%d ", k);
    printf("\teNodeB%d is seleced(%d).\n", i, associated_UE_Num);*/
    if(*(active + i) == 0){
      *(active + i) = bernoulli(*(binary_probability + i));
      if(*(active + i) == 1){/*
	for(j = 0; j < Ntp; j++){
	  if(*(ptrAdjacentBSs + i * Ntp + j) == 1)
	      printf("%d ", j);
	}
	printf("\n");
	for(j = 0; j < Ntp; j++){
	  if(*(associated_UE + j) == 0)
	    printf("%d ", j);
	}*/
	for(j = 0; j < Ntp; j++){
	  if(*(associated_UE + j) == 0 && *(ptrAdjacentBSs + i * Ntp + j) == 1){
	    *(associated_UE + j) = 1;
	    associated_UE_Num = associated_UE_Num + 1;
	    //printf("%d ", j);
	    
	  }
	}
	//printf("\neNodeB%d is active now with %d UEs.\n", i, associated_UE_Num);
      }
      else if(*(active + i) == -1){
	printf("operation mode error.\n");
	exit(1);
      }
    }
  }/*
  printf("%d\n", associated_UE_Num);
  for(i = 0; i < Tbs; i++)
    printf("%d\n", *(load_demand + i));*/
  free(associated_UE);
  free(binary_probability);
  free(load_demand);
  free(power);
  char name[80];
  char fileName[80];
  char str[10];
  sprintf(str, "%g", ABS);
  strcpy(name, "Static Activity Level Sleeping(");
  strcat(name, str);
  strcat(name, ")");
  strcpy(fileName, "StaticSleeping(");
  strcat(fileName, str);
  strcat(fileName, ").txt");
  FILE *PrintToFile;
  if((PrintToFile = fopen(fileName, "w")) == NULL)
    printf("\nerror! Fail to open %s!", fileName);
  else
    printf("\nOpen %s successfully!\n", fileName);
  fprintf(PrintToFile, "This is the input to CPLEX for %s.\n", name);
  cobjective(Tbs, PrintToFile, Ntp, Modulation, ptrBR, ptrDS);
  for(i = 0; i < Tbs - mNbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, 1 - ABS);
  for(i = Tbs - mNbs; i < Tbs; i++)
    fprintf(PrintToFile, "x%d - %g t = 0\n", i, ABS);
  ABSconstraint1(Tbs, ptrDS, PrintToFile, ptrminM, mNbs, Modulation, ptrDS, DSt, Ntp);
  cconstraint2(Tbs, PrintToFile, mNbs, Ntp, Modulation, ptrP, MP, mP);
  cconstraint3(Tbs, PrintToFile, nUser, Modulation);
  cconstraint4(Tbs, PrintToFile, Ntp, Modulation, BP);
  cconstraint5(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint6(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint7(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint8(Tbs, PrintToFile, Ntp, Modulation);
  cconstraint9(Tbs, PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, ptrP, mABP, mSBP);
  cconstraint10(Tbs, PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
  cconstraint11(Tbs, PrintToFile, nUser, Ntp, Modulation);
  for(i = 0; i < Tbs; i++)
    fprintf(PrintToFile, "b%d = %d\n", i, *(active + i));
  ABSbounds(Tbs, PrintToFile, Ntp, Modulation, mNbs);
  cspecifyTypes(Tbs, PrintToFile, Ntp, Modulation);
  fclose(PrintToFile);
  free(active);
}
