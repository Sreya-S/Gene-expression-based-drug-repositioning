#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

FILE *fp;

int a,b,i,j,k,l,m,n,k,s;

float xx;

char line[10000];

char strs[1000][1000];

double E[200][30000];

int affyN,GSMn,GN,GA[30000];
char G[30000][100],affys[30000][100];

/* non redundant gene set int and char */

int GGN;

char GG[30000][100];  

int CONTN,CONT[100],TREATN,TREAT[100];

double CONTAVE,TREATAVE;

double Xbar,Ybar,X2bar,Y2bar,XYbar,pears;

double sumsqrX,sumsqrY,Xbarsqr,Ybarsqr;

int R[100],RN;

double y[100];

int lineno,MMSE[100],BRAAK[100],AGE[100];

char str[100];

double pears,Z,Pcd[30000],Zcd[30000],Pmmse[30000],Zmmse[30000],Pbraak[30000],Zbraak[30000],Page[30000],Zage[30000];

int scores[30000],order[30000];

int MN,H[1000],sn,BR[1000],RPmmse[60000],RZmmse[60000];

double MMSEstattest();



int main()
{



if((fp=fopen("AtoG/AtoG-GPL96.txt","r"))==NULL){printf("no file\n");goto thend;}
affyN=0;
while(fgets(line,100,fp)!=NULL){
if(strlen(line)>0){
affyN=affyN+1; sscanf(line,"%s %d %s",&affys[affyN],&GA[affyN],&G[affyN]);
}}
fclose(fp);


lineno=-10000000;
fp=fopen("GSE/GSE1297_series_matrix.txt","r");
while(fgets(line,10000,fp)!=NULL){

lineno=lineno+1;

s=1;j=-1;for(i=0;i<=strlen(line);i++){	
j=j+1;
if((j>=0)&&(line[i]=='\t')){s=s+1;strs[s-1][j]='\0';j=-1;}
if(j>=0)strs[s][j]=line[i];
if(line[i]=='\n')strs[s][j]='\0';
}


j=-1;for(i=0;i<=strlen(strs[1]);i++){if(strs[1][i]!='\"'){j=j+1;strs[1][j]=strs[1][i];}}




/* reading characteristics */

if(strcmp("!Sample_characteristics_ch1",strs[1])==0){

for(k=1;k<=s;k++){
j=-1;for(i=0;i<=strlen(strs[k]);i++){if(strs[k][i]!='\"'){j=j+1;strs[k][j]=strs[k][i];}}	
}
	
sscanf(strs[2],"%s",&str);

if(strcmp(str,"mmse:")==0){
GSMn=0;
for(n=2;n<=s;n++){
sscanf(strs[n],"%s %d",&str,&m);
GSMn=GSMn+1;
MMSE[GSMn]=m;
}}


if(strcmp(str,"braak:")==0){
GSMn=0;
for(n=2;n<=s;n++){
sscanf(strs[n],"%s %d",&str,&b);
GSMn=GSMn+1;
BRAAK[GSMn]=b;
}}


if(strcmp(str,"age:")==0){
GSMn=0;
for(n=2;n<=s;n++){
sscanf(strs[n],"%s %d",&str,&a);
GSMn=GSMn+1;
AGE[GSMn]=a;
}}



}




if(strcmp(strs[1],"!series_matrix_table_end")==0){lineno=-1000000;}

if(strcmp(strs[1],"ID_REF")==0){lineno=0;}

if(lineno>0){

m=0;
for(i=1;i<=affyN;i++){if(strcmp(affys[i],strs[1])==0){m=i;goto POP;}}POP:;

if(m!=0){
GSMn=0;for(i=2;i<=s;i++){GSMn=GSMn+1;sscanf(strs[i],"%f",&xx);E[GSMn][m]=xx;}
}

}


}

fclose(fp);

for(n=1;n<=GSMn;n++)printf("%d %d %d %d\n",n,MMSE[n],AGE[n],BRAAK[n]);


MMSEstattest();goto thend;




fp=fopen("GSE1297everything.csv","w");

fprintf(fp,"Gene,cdPearson,cdZ,mPearson,mZ,bPearson,bZ,aPearson,aZ\n");

for(k=1;k<=affyN;k++){


/* MMSE vs gene expression*/

Xbar=0,Ybar=0,X2bar=0,Y2bar=0,XYbar=0;

for(n=1;n<=GSMn;n++){

Xbar=Xbar+E[n][k];
Ybar=Ybar+MMSE[n];
X2bar=X2bar+E[n][k]*E[n][k];
Y2bar=Y2bar+MMSE[n]*MMSE[n];
XYbar=XYbar+MMSE[n]*E[n][k];

}

pears=(XYbar-Xbar*Ybar/GSMn)/sqrt((Y2bar-Ybar*Ybar/GSMn)*(X2bar-Xbar*Xbar/GSMn));

Z=0.5*log((1+pears)/(1-pears))*sqrt(GSMn-3);

Pmmse[k]=pears;
Zmmse[k]=Z;

/* printf("%s,%f,%f\n",G[k],pears,Z); */

/*BRAAK vs gene expression */

Xbar=0,Ybar=0,X2bar=0,Y2bar=0,XYbar=0;

for(n=1;n<=GSMn;n++){

Xbar=Xbar+E[n][k];
Ybar=Ybar+BRAAK[n];
X2bar=X2bar+E[n][k]*E[n][k];
Y2bar=Y2bar+BRAAK[n]*BRAAK[n];
XYbar=XYbar+BRAAK[n]*E[n][k];
}

pears=(XYbar-Xbar*Ybar/GSMn)/sqrt((Y2bar-Ybar*Ybar/GSMn)*(X2bar-Xbar*Xbar/GSMn));

Z=0.5*log((1+pears)/(1-pears))*sqrt(GSMn-3);

Pbraak[k]=pears;
Zbraak[k]=Z;



/*AGE vs gene expression */

Xbar=0,Ybar=0,X2bar=0,Y2bar=0,XYbar=0;

for(n=1;n<=GSMn;n++){

Xbar=Xbar+E[n][k];
Ybar=Ybar+AGE[n];
X2bar=X2bar+E[n][k]*E[n][k];
Y2bar=Y2bar+AGE[n]*AGE[n];
XYbar=XYbar+AGE[n]*E[n][k];
}

pears=(XYbar-Xbar*Ybar/GSMn)/sqrt((Y2bar-Ybar*Ybar/GSMn)*(X2bar-Xbar*Xbar/GSMn));

Z=0.5*log((1+pears)/(1-pears))*sqrt(GSMn-3);

Page[k]=pears;
Zage[k]=Z;

fprintf(fp,"%s,%f,%f,%f,%f,%f,%f\n",G[k],Pmmse[k],Zmmse[k],Pbraak[k],Zbraak[k],Page[k],Zage[k]);

/* printf("%s,%f,%f\n",G[k],pears,Z); */
}

fclose(fp);





/* MMSE vs Gene Expression */


for(k=1;k<=affyN;k++)scores[k]=sqrt(Zmmse[k]*Zmmse[k]);


for(i=1;i<=affyN;i++)order[i]=i;
for(i=1;i<=affyN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(scores[k]>scores[order[j]]){order[j+1]=order[j];order[j]=k;}
}}


	

	
fp=fopen("GSE1297-HIPmmse.txt","w");

fprintf(fp,"Gene\tGene no.\tZscore\n");

GGN=0;
for(n=1;n<=affyN;n++){
i=order[n];

for(j=1;j<=GGN;j++){if(strcmp(G[i],GG[j])==0)goto LLLL;}

GGN=GGN+1;sprintf(GG[GGN],"%s",G[i]);

if(sqrt(Zmmse[i]*Zmmse[i])<2)goto LLLL;

GGN=GGN+1;

fprintf(fp,"%s\t%d\t%f\n",G[i],GA[i],Zmmse[i]);

LLLL:;
}


fclose(fp);



/* BRAAK vs Gene Expression */


for(k=1;k<=affyN;k++)scores[k]=sqrt(Zbraak[k]*Zbraak[k]);


for(i=1;i<=affyN;i++)order[i]=i;
for(i=1;i<=affyN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(scores[k]>scores[order[j]]){order[j+1]=order[j];order[j]=k;}
}}

/* for(i=1;i<=10;i++)printf("%s,%f\n",G[order[i]],Zbraak[order[i]]); */

fp=fopen("GSE1297-HIPbraak.txt","w");

fprintf(fp,"Gene\tGene no.\tZscore\n");

GGN=0;
for(n=1;n<=affyN;n++){
i=order[n];

for(j=1;j<=GGN;j++){if(strcmp(G[i],GG[j])==0)goto MMMM;}

GGN=GGN+1;sprintf(GG[GGN],"%s",G[i]);

if(sqrt(Zbraak[i]*Zbraak[i])<2)goto MMMM;

GGN=GGN+1;

fprintf(fp,"%s\t%d\t%f\n",G[i],GA[i],Zbraak[i]);

MMMM:;
}


fclose(fp);



/* AGE vs Gene Expression */


for(k=1;k<=affyN;k++)scores[k]=sqrt(Zage[k]*Zage[k]);


for(i=1;i<=affyN;i++)order[i]=i;
for(i=1;i<=affyN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(scores[k]>scores[order[j]]){order[j+1]=order[j];order[j]=k;}
}}

/* for(i=1;i<=10;i++)printf("%s,%f\n",G[order[i]],Zage[order[i]]); */

fp=fopen("GSE1297-HIPage.txt","w");

fprintf(fp,"Gene\tGene no.\tZscore\n");

GGN=0;
for(n=1;n<=affyN;n++){
i=order[n];

for(j=1;j<=GGN;j++){if(strcmp(G[i],GG[j])==0)goto NNNN;}

GGN=GGN+1;sprintf(GG[GGN],"%s",G[i]);

if(sqrt(Zage[i]*Zage[i])<2)goto NNNN;

GGN=GGN+1;

fprintf(fp,"%s\t%d\t%f\n",G[i],GA[i],Zage[i]);

NNNN:;
}


fclose(fp);






thend:;







return 0; 
}

double MMSEstattest()
{
	
int Z2n,Z2nR,count,randN;

double threshold=2;

/**Random data**/

srand((unsigned int)time(NULL));

/* restrict to hippocampus with mmse readings */
MN=0;
for(i=1;i<=GSMn;i++){
if(MMSE[i]!=-999){
MN=MN+1;H[MN]=i;
}}
for(i=1;i<=MN;i++)printf("%d %d %d\n",i,H[i],MMSE[H[i]]);
/**************************************************/

Z2n=0;
for(k=1;k<=affyN;k++){
	
Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;n=0;

for(j=1;j<=MN;j++){
i=H[j];
n=n+1;
Xbar=Xbar+E[i][k];
Ybar=Ybar+MMSE[i];
X2bar=X2bar+E[i][k]*E[i][k];
Y2bar=Y2bar+MMSE[i]*MMSE[i];
XYbar=XYbar+E[i][k]*MMSE[i];
}

pears=(XYbar-Xbar*Ybar/n)/sqrt((Y2bar-Ybar*Ybar/n)*(X2bar-Xbar*Xbar/n));

Z=0.5*log((1+pears)/(1-pears))*sqrt(n-3);

if(sqrt(Z*Z)>threshold)Z2n=Z2n+1;

}


for(randN=1;randN<=1000;randN++){


/* random permutation of the set order[i] is random perm of i ********************/
for(i=1;i<=MN;i++)R[i]=rand()%1000+1;

for(i=1;i<=MN;i++)order[i]=i;
for(i=1;i<=MN-1;i++){k=order[i+1];
for(j=i;j>=1;j--){
if(R[k]>R[order[j]]){order[j+1]=order[j];order[j]=k;}
}}
/******************************************************/

Z2nR=0;
for(k=1;k<=affyN;k++){




/* output Z scores */


Xbar=0;Ybar=0;X2bar=0;Y2bar=0;XYbar=0;n=0;

for(j=1;j<=MN;j++){
i=H[j];m=H[order[j]];
n=n+1;
Xbar=Xbar+E[i][k];
Ybar=Ybar+MMSE[m];
X2bar=X2bar+E[i][k]*E[i][k];
Y2bar=Y2bar+MMSE[m]*MMSE[m];
XYbar=XYbar+E[i][k]*MMSE[m];
}


pears=(XYbar-Xbar*Ybar/n)/sqrt((Y2bar-Ybar*Ybar/n)*(X2bar-Xbar*Xbar/n));

Z=0.5*log((1+pears)/(1-pears))*sqrt(n-3);

if(sqrt(Z*Z)>threshold)Z2nR=Z2nR+1;


}	

if(Z2nR>Z2n)count=count+1;

}

printf("probability = %f",(double)count/1000);



return 0;
}

