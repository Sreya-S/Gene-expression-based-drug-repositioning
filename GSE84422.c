#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

char fileo[100];

int a,b,i,j,k,l,m,n,k,s,h;

char line[10000];

char strs[1500][1000];

double E[1500][60000]; /****Expression has sample numbers and gene values***/
int affyN,GSMn,GN,GA[60000];
char G[60000][100],affys[60000][100];

int MMSE[1500],AGE[1500];

double BRAAK[1500],APOE[1500],CERAD[1500],CDR[1500],PLAQ[1500],TANGLES[1500];

int yes;

/* non redundant gene set int and char */
char GG[60000][100];

int lineno;

char str[1500],str1[1500];
int GGN;  

int sn;

int SEX[300],BR[1500];

double Xbar,Ybar,X2bar,Y2bar,XYbar,pears;

double sumsqrX,sumsqrY,Xbarsqr,Ybarsqr;

double y[300];

double pears,Z,PP[60000],ZZ[60000];

int scores[60000],order[60000];

int R[1500];

double PearsCalculation();

double Y[1500];
int II[1500];

int N;

double RANDOMtest();
double NonRedundant();


double F[1000][60000];

int BRn=19;
int br;
char BRAIN[20][4]={"XXX","FP","OVC","ITG","MTG","STG","PCC","AC","PHG","TP","PG","IFG","DPC","SPL","PC","CN","HIP","PUT","AMG","NA"};

char GSE[1000],GSEsamples[1000],GPL[100],AtoG[100],OUTPUT[100];

int main()
{



	
sprintf(GPL,"GPL96");
sprintf(AtoG,"AtoG/AtoG-%s.txt",GPL);
sprintf(GSE,"GSE/GSE84422-%s_series_matrix.txt",GPL);
sprintf(GSEsamples,"GSE/GSE84422-%s-SAMPLES.txt",GPL);


fp=fopen(GSEsamples,"r");
fgets(line,10000,fp);
while(fgets(line,10000,fp)!=NULL){
	
s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';	
}



sscanf(strs[14],"%d",&sn);

for(i=1;i<=BRn;i++){if(strcmp(strs[13],BRAIN[i])==0){BR[sn]=i;goto QQQQ;}}QQQQ:;

sscanf(strs[2],"%d",&AGE[sn]);

if(strcmp(strs[3],"male")==0)SEX[sn]=1;
if(strcmp(strs[3],"female")==0)SEX[sn]=2;

sscanf(strs[8],"%lf",&BRAAK[sn]);
sscanf(strs[11],"%lf",&CERAD[sn]);
sscanf(strs[7],"%lf",&CDR[sn]);
sscanf(strs[10],"%lf",&PLAQ[sn]);
sscanf(strs[12], "%lf",&TANGLES[sn]);

}
fclose(fp);	

if((fp=fopen(AtoG,"r"))==NULL)

affyN=0;
while(fgets(line,1000,fp)!=NULL){
if(strlen(line)>0){
affyN=affyN+1; sscanf(line,"%s %d %s",&affys[affyN],&GA[affyN],&G[affyN]);
}}
fclose(fp);

/*printf("affyN %d\n",affyN);*/

lineno=-10000000;
fp=fopen(GSE,"r");
while(fgets(line,10000,fp)!=NULL){

lineno=lineno+1;

s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';
}

j=-1;for(i=0;i<=strlen(strs[1]);i++){if(strs[1][i]!='\"'){j=j+1;strs[1][j]=strs[1][i];}}


if(strcmp(strs[1],"ID_REF")==0){lineno=0;}

if(lineno>0){	


m=0;
for(i=1;i<=affyN;i++){if(strcmp(affys[i],strs[1])==0){m=i;goto POP;}}POP:;


if(m!=0){
GSMn=0;for(i=2;i<=s;i++){GSMn=GSMn+1;sscanf(strs[i],"%lf",&E[GSMn][m]);}
}

}

}

fclose(fp);  

/*
fp=fopen("test.txt","w");
for(i=1;i<=affyN;i++){
fprintf(fp,"%s\t",affys[i]);
for(j=1;j<=GSMn;j++)fprintf(fp,"%f\t",E[j][i]);
fprintf(fp,"\n");
}
fclose(fp);
*/



for(br=1;br<=BRn;br++){
	

N=0;
for(i=1;i<=sn;i++){
if(BR[i]==br){
N=N+1;Y[N]=CERAD[i];II[N]=i;
}
}

if(N==0)goto AAAA;

RANDOMtest();

sprintf(OUTPUT,"PROFILES/GSE8442296-%svsCERAD.txt",BRAIN[br]);

fp=fopen(OUTPUT,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);


/* Non-redundant set */

sprintf(OUTPUT,"PROFILES/GSE8442296-%svsCERAD-nonredundant.txt",BRAIN[br]);

fp=fopen(OUTPUT,"w");

NonRedundant();

fclose(fp);

AAAA:;
}



for(br=1;br<=BRn;br++){
	

N=0;
for(i=1;i<=sn;i++){
if(BR[i]==br){
N=N+1;Y[N]=BRAAK[i];II[N]=i;
}
}

if(N==0)goto BBBB;

RANDOMtest();

sprintf(OUTPUT,"PROFILES/GSE8442296-%svsBRAAK.txt",BRAIN[br]);

fp=fopen(OUTPUT,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);


/* Non-redundant set */

sprintf(OUTPUT,"PROFILES/GSE8442296-%svsBRAAK-nonredundant.txt",BRAIN[br]);

fp=fopen(OUTPUT,"w");

NonRedundant();

fclose(fp);

BBBB:;
}











return 0;
}

double RANDOMtest()

{
	
double Zcut=2;
int N1,Nr,randN,M;

srand((unsigned int)time(NULL));

for(k=1;k<=affyN;k++){

pears=0;Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;

for(i=1;i<=N;i++){
Xbar=Xbar+E[II[i]][k];
Ybar=Ybar+Y[i];
X2bar=X2bar+E[II[i]][k]*E[II[i]][k];
Y2bar=Y2bar+Y[i]*Y[i];
XYbar=XYbar+E[II[i]][k]*Y[i];
}


if((Y2bar-Ybar*Ybar/N)*(X2bar-Xbar*Xbar/N)>=0){
pears=(XYbar-Xbar*Ybar/N)/sqrt((Y2bar-Ybar*Ybar/N)*(X2bar-Xbar*Xbar/N));
}

Z=0.5*log((1+pears)/(1-pears))*sqrt(N-3);

ZZ[k]=Z;
PP[k]=pears;

if(sqrt(Z*Z)>Zcut)N1=N1+1;

}

goto SKIPRAND;


M=0;
for(randN=1;randN<=1000;randN++){


/* random permutation of the set order[i] is random perm of i ********************/
for(i=1;i<=N;i++)R[i]=rand()%1000+1;

for(i=1;i<=N;i++)order[i]=i;
for(i=1;i<=N-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(R[k]>R[order[j]]){order[j+1]=order[j];order[j]=k;}
}}



Nr=0;
for(k=1;k<=affyN;k++){

pears=0;Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;

for(i=1;i<=N;i++){
Xbar=Xbar+E[II[i]][k];
Ybar=Ybar+Y[R[i]];
X2bar=X2bar+E[II[i]][k]*E[II[i]][k];
Y2bar=Y2bar+Y[R[i]]*Y[R[i]];
XYbar=XYbar+E[II[i]][k]*Y[R[i]];
}


if((Y2bar-Ybar*Ybar/N)*(X2bar-Xbar*Xbar/N)>=0){
pears=(XYbar-Xbar*Ybar/N)/sqrt((Y2bar-Ybar*Ybar/N)*(X2bar-Xbar*Xbar/N));
}

Z=0.5*log((1+pears)/(1-pears))*sqrt(N-3);

if(sqrt(Z*Z)>Zcut)Nr=Nr+1;

}
printf("%d\n",Nr);


if(Nr>N1)M=M+1;

}

printf("%s %d %d %d %f\n",BRAIN[br],N1,M,randN,(double)M/(double)randN);


SKIPRAND:;

return 0;
}




double NonRedundant()

{

for(k=1;k<=affyN;k++)scores[k]=sqrt(ZZ[k]*ZZ[k]);


for(i=1;i<=affyN;i++)order[i]=i;
for(i=1;i<=affyN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(scores[k]>scores[order[j]]){order[j+1]=order[j];order[j]=k;}
}}

/* for(i=1;i<=10;i++)printf("%s,%f\n",G[order[i]],ZZ[order[i]]); */
	

fprintf(fp,"Gene\tGene no.\tZscore\n");

GGN=0;
for(n=1;n<=affyN;n++){
i=order[n];

for(j=1;j<=GGN;j++){if(strcmp(G[i],GG[j])==0)goto SKIP;}

GGN=GGN+1;sprintf(GG[GGN],"%s",G[i]);

if(sqrt(ZZ[i]*ZZ[i])<2)goto SKIP;

fprintf(fp,"%s\t%d\t%f\n",G[i],GA[i],ZZ[i]);

SKIP:;

}


return 0;

}
























