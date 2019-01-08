#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;
FILE *fp1;

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

double NN,Xbar,Ybar,X2bar,Y2bar,XYbar,pears;

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
double ReadingFiles();


double F[1000][60000];

int BRn=19;
int br;
char BRAIN[20][4]={"XXX","FP","OVC","ITG","MTG","STG","PCC","AC","PHG","TP","PG","IFG","DPC","SPL","PC","CN","HIP","PUT","AMG","NA"};

char GSE[1000],GSEsamples[1000],GPL[100],AtoG[100],OUTPUT[100];

char FILES[200][100],NAMES[200][100];

int FILESn,fn;

int Qp,Qm;

int main()
{
	
	
GN=0;


FILESn=0;

fp=fopen("GSE/SETS-GSE1297BRAAKvREST.txt","r");
while(fgets(line,1000,fp)!=NULL){
sscanf(line,"%d %s %s",&i,&strs[1],&strs[2]);
FILESn=FILESn+1;
sprintf(NAMES[FILESn],"%s",strs[1]);
sprintf(FILES[FILESn],"%s",strs[2]);printf("%d %s %s\n",FILESn,NAMES[FILESn],FILES[FILESn]);
}
fclose(fp);


for(fn=1;fn<=FILESn;fn++){
ReadingFiles();
m=0;for(k=1;k<=GN;k++){if(F[fn][k]!=0)m=m+1;}
printf("%s %d\n",NAMES[fn],m);
}

fp=fopen("GSE1297BRAAKvREST-TABLE1.csv","w");

i=1;

n=0;
for(k=1;k<=GN;k++){if(F[i][k]!=0)n=n+1;}
printf("number of genes %d\n",n);
	
for(j=2;j<=FILESn;j++){


Z=0;pears=0;NN=0;Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;

for(k=1;k<=GN;k++){
if(F[i][k]*F[j][k]!=0){
NN=NN+1;
Xbar=Xbar+F[i][k];
Ybar=Ybar+F[j][k];
X2bar=X2bar+F[i][k]*F[i][k];
Y2bar=Y2bar+F[j][k]*F[j][k];
XYbar=XYbar+F[i][k]*F[j][k];
}}



if(NN<4)goto SKIP;

if((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN)>=0){
pears=(XYbar-Xbar*Ybar/NN)/sqrt((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN));
}

if(pears*pears==1)pears=0.99*pears;

Z=0.5*log((1+pears)/(1-pears))*sqrt(NN-3);

SKIP:;

fprintf(fp,"%s,%2.2f,%2.2f\n",NAMES[j],Z,pears);

}

fclose(fp);

fp=fopen("theBRAAKprofile.csv","w");

for(i=1;i<=GN;i++){

Qp=0;Qm=0;
for(j=1;j<=FILESn;j++){if(F[j][i]>0)Qp=Qp+1;if(F[j][i]<0)Qm=Qm+1;}

if(Qp+Qm>0)fprintf(fp,"%s,%d,%d\n",G[i],Qp,Qm);
}

fclose(fp);


goto thend;




fp1=fopen("1297mmseGRAPH.txt","w");

fp=fopen("1297mmseMATRIX.csv","w");

fprintf(fp,"NAMES,");for(i=1;i<=FILESn;i++)fprintf(fp,"%s,",NAMES[i]);fprintf(fp,"\n");

for(i=1;i<=FILESn;i++){
	
fprintf(fp,"%s,",NAMES[i]);
	
for(j=1;j<=FILESn;j++){


Z=0;pears=0;NN=0;Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;

for(k=1;k<=GN;k++){
if(F[i][k]*F[j][k]!=0){
NN=NN+1;
Xbar=Xbar+F[i][k];
Ybar=Ybar+F[j][k];
X2bar=X2bar+F[i][k]*F[i][k];
Y2bar=Y2bar+F[j][k]*F[j][k];
XYbar=XYbar+F[i][k]*F[j][k];
}}

if(NN<4)goto SKIP1;

if((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN)>=0){
pears=(XYbar-Xbar*Ybar/NN)/sqrt((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN));
}

if(pears*pears==1)pears=0.99*pears;

Z=0.5*log((1+pears)/(1-pears))*sqrt(NN-3);

SKIP1:;

fprintf(fp,"%2.2f,",pears);

if((i>1)&&(j>1)){
if(i!=j)fprintf(fp1,"%d\t%d\t%f\n",i,j,Z);
else fprintf(fp1,"%d\t%d\t%f\n",i,j,0);
}

}
fprintf(fp,"\n");
}

fclose(fp);

fclose(fp1);


fp=fopen("1297mmse.txt","w");
/*
fprintf(fp,"set xtics (");
for(i=1;i<=FILESn-1;i++)fprintf(fp,"\"%s\",",NAMES[i]);
fprintf(fp,"\"%s\")\n",NAMES[i]);
fprintf(fp,"set ytics rotate by 90\n");
fprintf(fp,"set ytics (");
for(i=1;i<=FILESn-1;i++)fprintf(fp,"\"%s\",",NAMES[i]);
fprintf(fp,"\"%s\")\n",NAMES[i]);
*/

/*fprintf(fp,"set xtics 1");*/

fprintf(fp,"plot \"1297mmseGRAPH.txt\" using 1:2 with lines notitle\n");
fclose(fp);

system("\"\"C:/Program Files (x86)/gnuplot/bin/gnuplot\" -persist \"1297mmse.txt\"\"");



thend:;

return 0;
	
}



double ReadingFiles()
{


for(k=1;k<=50000;k++)F[fn][k]=0;	

if((fp=fopen(FILES[fn],"r"))!=NULL){

fgets(line,1000,fp);
while(fgets(line,1000,fp)!=NULL){

sscanf(line,"%s %d %lf",&str,&k,&Z); if(k>GN)GN=k;

F[fn][k]=Z;sprintf(G[k],"%s",str);



}

fclose(fp);	
n=0;
for(k=1;k<=50000;k++){if(F[fn][k]!=0)n=n+1;}
printf("%s %d\n",FILES[fn],n);
}

return 0;		
}
























