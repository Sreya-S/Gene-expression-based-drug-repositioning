#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

int a,b,g,i,j,k,l,m,n,k,s;

float xx;

char line[200000];

char strs[1000][1000];

double E[1000][30000];

int affyN,GSMn,GN,GA[30000];
char G[30000][100],affys[30000][100];

/* non redundant gene set int and char */

int GGN;

int GG[30000];  


double XYn,Xbar,Ybar,X2bar,Y2bar,XYbar,pears,Z;


int R[100],RN;

double y[100];

int lineno,GROUP[300];

char str[100];

double PP[50000],ZZ[50000];

int scores[30000],order[30000];

int sn,TYPE[1000],RPmmse[60000],RZmmse[60000];

double MMSEstattest();

int CMPn;
char CMP[2000][100];
int CMPi[2000];

int TREAT[2000];

int cmpN;

int maxi;

int main()
{

CMPn=0;
GSMn=0;
fp=fopen("E:/3rd year/Sreya/PROGRAMS/GSE/GSE85871-SAMPLES.txt","r");

while(fgets(line,10000,fp)!=NULL){
	
sscanf(line,"%s %s",&strs[1],&strs[2]);	

GSMn=GSMn+1;
	
for(i=1;i<=CMPn;i++){if(strcmp(strs[2],CMP[i])==0){CMPi[GSMn]=i;goto QQQQ;}}
CMPn=CMPn+1;sprintf(CMP[CMPn],"%s",strs[2]);CMPi[GSMn]=CMPn;
QQQQ:;
	
	
	
}

fclose(fp);

for(i=1;i<=GSMn;i++)printf("%d %d %s\n",i,CMPi[i],CMP[CMPi[i]]);

/* DMSO #31 is our control */


	
GN=0;	
fp=fopen("E:/3rd year/Sreya/PROGRAMS/AtoG/AtoG-GPL571.txt","r");
affyN=0;
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
affyN=affyN+1; 
sscanf(line,"%s %d %s",&affys[affyN],&GA[affyN],&G[affyN]);if(GA[affyN]>GN)GN=GA[affyN];
}}
fclose(fp);

printf("affyN %d\n",affyN);



maxi=0;
lineno=-10000000;
fp=fopen("E:/3rd year/Sreya/PROGRAMS/GSE/GSE85871_series_matrix.txt","r");
while(fgets(line,100000,fp)!=NULL){
	
/*if(strlen(line)>maxi){maxi=strlen(line);printf("%d\n",strlen(line));} */

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




fp=fopen("DATABASE1.txt","w");

	
for(cmpN=1;cmpN<=CMPn;cmpN++){
if(cmpN!=31){
	
printf("%s\n",CMP[cmpN]);
	
for(n=1;n<=GSMn;n++){
TREAT[n]=0;
if(CMPi[n]==cmpN)TREAT[n]=2;
if(CMPi[n]==31)TREAT[n]=1;
}

m=0;
for(k=1;k<=affyN;k++){
	
XYn=0;Z=0;pears=0;XYn=0;Xbar=0,Ybar=0,X2bar=0,Y2bar=0,XYbar=0;

for(n=1;n<=GSMn;n++){
if(TREAT[n]>0){
XYn=XYn+1;

Xbar=Xbar+E[n][k];
Ybar=Ybar+TREAT[n];
X2bar=X2bar+E[n][k]*E[n][k];
Y2bar=Y2bar+TREAT[n]*TREAT[n];
XYbar=XYbar+TREAT[n]*E[n][k];

}}

if((Y2bar-Ybar*Ybar/XYn)*(X2bar-Xbar*Xbar/XYn)>0){
pears=(XYbar-Xbar*Ybar/XYn)/sqrt((Y2bar-Ybar*Ybar/XYn)*(X2bar-Xbar*Xbar/XYn));

if(XYn>3)Z=0.5*log((1+pears)/(1-pears))*sqrt(XYn-3);
}


if(sqrt(Z*Z)<=2)Z=0;
else m=m+1;

ZZ[k]=Z;


}

printf("%d\n",m);


for(k=1;k<=affyN;k++)scores[k]=sqrt(ZZ[k]*ZZ[k]);


for(i=1;i<=affyN;i++)order[i]=i;
for(i=1;i<=affyN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(scores[k]>scores[order[j]]){order[j+1]=order[j];order[j]=k;}
}}
 
	
GGN=0;
for(n=1;n<=affyN;n++){
i=order[n];

for(j=1;j<=GGN;j++){if(GA[i]==GG[j])goto SKIP;}

if(ZZ[i]==0)goto SKIP;

GGN=GGN+1;GG[GGN]=GA[i];

PP[GGN]=ZZ[i];

SKIP:;

}

fprintf(fp,">%s %d\n",CMP[cmpN],GGN);
for(i=1;i<=GGN;i++)fprintf(fp,"%d %2.2f ",GG[i],PP[i]);
fprintf(fp,"\n");

}}
	
fclose(fp);

	

thend:;	
	
	
	
	
return 0;
}
