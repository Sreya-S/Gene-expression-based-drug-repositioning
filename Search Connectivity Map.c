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



double NN,Xbar,Ybar,X2bar,Y2bar,XYbar,pears,Z;
double NNI,XbarI,YbarI,X2barI,Y2barI,XYbarI,pearsI,ZI;

double QI,CI;


double PP[60000],ZZ[60000];

double PPI[60000],ZZI[60000];

int scores[60000],order[60000];

int CMPn;
char CMP[10000][100];

int N;

double x;


int FILESn,fn;

double QUERY[30000];

double F[1500][30000];

int Qp,Qm;

int main()
{
	
GN=29999;
	
/* read in the query */

for(i=1;i<=GN;i++)QUERY[i]=0;

fp=fopen("PROFILES/GSE8442296-HIPvsBRAAK96-nonredundant.txt","r");
fgets(line,100,fp);
while(fgets(line,100,fp)!=NULL){

sscanf(line,"%s %d %lf",&str,&k,&x);

QUERY[k]=x;sprintf(G[k],"%s",str);
}

fclose(fp);

Qp=0;Qm=0;
for(i=1;i<=GN;i++){if(QUERY[i]<0)Qm=Qm+1;if(QUERY[i]>0)Qp=Qp+1;}

printf("%d %d\n",Qp,Qm);

/* read database */

fp=fopen("ConnectivityMap/HGNC-CMAP.txt","r");
CMPn=0;
for(fn=1;fn<=100000;fn++){

fscanf(fp,"%s",&str);if(feof(fp))goto WWWW;	
fscanf(fp,"%d",&k);if(feof(fp))goto WWWW;	

CMPn=CMPn+1;
sprintf(CMP[CMPn],"%s",str);

for(i=1;i<=GN;i++)F[CMPn][i]=0;

Z=0;pears=0;NN=0;Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;

ZI=0;pearsI=0;NNI=0;XbarI=0;YbarI=0;X2barI=0;Y2barI=0;XYbarI=0;

for(i=1;i<=k;i++){fscanf(fp,"%d",&j);fscanf(fp,"%lf",&x);F[CMPn][j]=x;

if(QUERY[j]!=0){
NN=NN+1;
Xbar=Xbar+QUERY[j];
Ybar=Ybar+x;
X2bar=X2bar+QUERY[j]*QUERY[j];
Y2bar=Y2bar+x*x;
XYbar=XYbar+QUERY[j]*x;

QI=QUERY[i]/sqrt(QUERY[i]*QUERY[i]);
CI=x/sqrt(x*x);
NNI=NNI+1;
XbarI=XbarI+QI;
YbarI=YbarI+CI;
X2barI=X2barI+QI*QI;
Y2barI=Y2barI+CI*CI;
XYbarI=XYbarI+QI*CI;

}}

if(NN<4)goto SKIP;

if((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN)>=0){
pears=(XYbar-Xbar*Ybar/NN)/sqrt((Y2bar-Ybar*Ybar/NN)*(X2bar-Xbar*Xbar/NN));
}

Z=0.5*log((1+pears)/(1-pears))*sqrt(NN-3);


if((Y2barI-YbarI*YbarI/NNI)*(X2barI-XbarI*XbarI/NNI)>=0){
pearsI=(XYbarI-XbarI*YbarI/NNI)/sqrt((Y2barI-YbarI*YbarI/NNI)*(X2barI-XbarI*XbarI/NNI));
}

ZI=0.5*log((1+pearsI)/(1-pearsI))*sqrt(NNI-3);


SKIP:;

ZZ[CMPn]=Z;

ZZI[CMPn]=ZI;

}

WWWW:;
fclose(fp);


for(i=1;i<=CMPn;i++)order[i]=i;
for(i=1;i<=CMPn-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(ZZ[k]<ZZ[order[j]]){order[j+1]=order[j];order[j]=k;}
}}

for(i=1;i<=20;i++){j=order[i];printf("%s %f %d %f\n",CMP[j],ZZ[j],i,ZZI[j]);}

for(n=1;n<=20;n++){
k=order[n];
sprintf(fileo,"cMAP-84422/cMAP-84422x-hipGENESreversedBY-%s.csv",CMP[k]);
printf("%s\n",fileo);
fp=fopen(fileo,"w");
for(i=1;i<=GN;i++){
if(QUERY[i]*F[k][i]<0){
fprintf(fp,"%s,%f,%f\n",G[i],QUERY[i],F[k][i]);
}
}}
fclose(fp);
	
thend:;

return 0;
}
