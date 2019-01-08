#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

char fileo[100];

int a,b,i,j,k,l,m,n,k,s,h;

char line[10000];

char strs[1000][1000];

double E[1000][60000]; /****Expression has sample numbers and gene values***/
int affyN,GSMn,GN,GA[60000];
char G[60000][1000],affys[60000][1000];

int MMSE[300],AGE[300];

double BRAAK[300],APOE[300];

int yes;

/* non redundant gene set int and char */
int GG[60000][1000];

int lineno;

char str[1000],str1[1000];
int GGN;  

int sn;

int SEX[300],BR[1000];

double Xbar,Ybar,X2bar,Y2bar,XYbar,pears;

double sumsqrX,sumsqrY,Xbarsqr,Ybarsqr;

double y[300];

double pears,Z,PP[60000],ZZ[60000];

int scores[60000],order[60000];

int R[1000];

double PearsCalculation();

double Y[1000];
int II[1000];

int N;

double RANDOMtest();
double NonRedundant();


int main()
{
	

	
fp=fopen("GSE/GSE48350_descriptors-new.txt","r");
fgets(line,10000,fp);
while(fgets(line,10000,fp)!=NULL){


s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';	
}

sscanf(strs[8],"%d",&sn);

if(strcmp(strs[1],"PCG")==0)BR[sn]=1;
if(strcmp(strs[1],"SFG")==0)BR[sn]=2;
if(strcmp(strs[1],"HIP")==0)BR[sn]=3;
if(strcmp(strs[1],"EC")==0)BR[sn]=4;  

if(strcmp(strs[3],"male")==0)SEX[sn]=1;
if(strcmp(strs[3],"female")==0)SEX[sn]=2;

sscanf(strs[4],"%d",&AGE[sn]);

sscanf(strs[5],"%lf",&BRAAK[sn]);
sscanf(strs[6],"%lf",&APOE[sn]);
sscanf(strs[7],"%d",&MMSE[sn]);


}
fclose(fp);	



if((fp=fopen("AtoG/AtoG-GPL570.txt","r"))==NULL)

affyN=0;
while(fgets(line,1000,fp)!=NULL){
if(strlen(line)>0){
affyN=affyN+1; sscanf(line,"%s %d %s",&affys[affyN],&GA[affyN],&G[affyN]);
}}
fclose(fp);


lineno=-10000000;
fp=fopen("GSE/GSE48350_series_matrix.txt","r");
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



/* MMSE & PCG */

sprintf(fileo,"MMSE&PCG.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==1)&&(MMSE[i]!=-999)){
N=N+1;Y[N]=MMSE[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* MMSE & PCG Ordered (Non-redundant) */

	
fp=fopen("GSE48350_PCG_MMSE.txt","w");

NonRedundant();

fclose(fp);


/* BRAAK & PCG */

sprintf(fileo,"BRAAK&PCG.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==1)&&(BRAAK[i]!=-999)){
N=N+1;Y[N]=BRAAK[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* BRAAK & PCG Ordered (Non-redundant) */

	
fp=fopen("GSE48350_PCG_BRAAK.txt","w");

NonRedundant();

fclose(fp);


/* MMSE & SFG */

sprintf(fileo,"MMSE&SFG.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==2)&&(MMSE[i]!=-999)){
N=N+1;Y[N]=MMSE[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* MMSE & SFG Ordered (Non-redundant) */

	
fp=fopen("GSE48350_SFG_MMSE.txt","w");

NonRedundant();

fclose(fp);

/* BRAAK & SFG */

sprintf(fileo,"BRAAK&SFG.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==2)&&(BRAAK[i]!=-999)){
N=N+1;Y[N]=BRAAK[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* BRAAK & SFG Ordered (Non-redundant) */

	
fp=fopen("GSE48350_SFG_BRAAK.txt","w");

NonRedundant();

fclose(fp);


/* MMSE & HIP */

sprintf(fileo,"MMSE&HIP.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==3)&&(MMSE[i]!=-999)){
N=N+1;Y[N]=MMSE[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* MMSE & HIP Ordered (Non-redundant) */

	
fp=fopen("GSE48350_HIP_MMSE.txt","w");

NonRedundant();

fclose(fp);

/* BRAAK & HIP */

sprintf(fileo,"BRAAK&HIP.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==3)&&(BRAAK[i]!=-999)){
N=N+1;Y[N]=BRAAK[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* BRAAK & HIP Ordered (Non-redundant) */

	
fp=fopen("GSE48350_HIP_BRAAK.txt","w");

NonRedundant();

fclose(fp);

/* MMSE & EC */

sprintf(fileo,"MMSE&EC.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==4)&&(MMSE[i]!=-999)){
N=N+1;Y[N]=MMSE[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);


/* MMSE & EC Ordered (Non-redundant) */

	
fp=fopen("GSE48350_EC_MMSE.txt","w");

NonRedundant();

fclose(fp);

/* BRAAK & EC */

sprintf(fileo,"BRAAK&EC.txt");

N=0;
for(i=1;i<=sn;i++){
if((BR[i]==4)&&(BRAAK[i]!=-999)){
N=N+1;Y[N]=BRAAK[i];II[N]=i;
}
}


RANDOMtest();


fp=fopen(fileo,"w");
for(k=1;k<=affyN;k++){
fprintf(fp,"%s\tGene= %s\tZ= %f\tPears= %f\n",affys[k],G[k],ZZ[k],PP[k]); /*genes not affys*/
}
fclose(fp);

/* BRAAK & EC Ordered (Non-redundant) */

	
fp=fopen("GSE48350_EC_BRAAK.txt","w");

NonRedundant();

fclose(fp);

	
thend:;

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




M=0;
for(randN=1;randN<=1000;randN++){


/* random permutation of the set order[i] is random perm of i ********************/
for(i=1;i<=N;i++)R[i]=rand()%1000+1;

for(i=1;i<=N;i++)order[i]=i;
for(i=1;i<=N-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(R[k]>R[order[j]]){order[j+1]=order[j];order[j]=k;}
}}
/******************************************************/


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

if(Nr>N1)M=M+1;

}

printf("%d %d %d %f\n",N1,M,randN,(double)M/(double)randN);

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

GGN=GGN+1;

fprintf(fp,"%s\t%d\t%f\n",G[i],GA[i],ZZ[i]);

SKIP:;

}


fclose(fp);
}







double PearsCalculation()

{

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

}

return 0;
}




