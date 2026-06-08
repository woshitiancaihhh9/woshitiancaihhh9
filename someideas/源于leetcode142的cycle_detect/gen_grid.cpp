#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
inline int Npos(int x,int mu,int lam){return x<mu?x:mu+((x-mu)%lam);}
inline bool eqp(int x,int y,int mu,int lam){ if(x<mu&&y<mu) return x==y; return Npos(x,mu,lam)==Npos(y,mu,lam);}

struct RA{ ll det, ent; int o; int mode; };
ll locate_cost(int oldc,int l,int mu,int lam){ ll en=l; int a=oldc,b=oldc+l; while(true){ if(eqp(a,b,mu,lam)) return en; en+=2; a++; b++; } }
RA A(int mu,int lam){ ll nextc=0; int oldc=0,c=0,step=1,skipEnd=0; bool skip=false;
  while(true){ int p; if(skip){p=skipEnd; skip=false;} else { p=c; for(int s=1;s<=step;s++){ nextc++; p++; if(eqp(p,c,mu,lam)){ ll e=locate_cost(oldc,s,mu,lam); return {nextc,e,oldc,0};} } }
    int newc=p; for(int t=1;t<=2*step;t++){ nextc++; p++; if(eqp(p,c,mu,lam)){ ll e=locate_cost(oldc,step+t,mu,lam); return {nextc,e,oldc,0};} if(eqp(p,newc,mu,lam)){ oldc=c; ll e=locate_cost(oldc,t,mu,lam); return {nextc,e,oldc,1};} }
    c=newc; step*=2; skipEnd=p; skip=true; } }

struct C{ll n,c;};
C Bsim(int mu,int lam){ ll nextc=1,cmp=0; int power=1,l=1,tort=0,hare=1; ll det=0; while(true){cmp++; if(eqp(tort,hare,mu,lam)) break; if(power==l){tort=hare; power*=2; l=0;} hare++; nextc++; l++;} det=nextc; tort=0; hare=l; nextc+=l; while(true){cmp++; if(eqp(tort,hare,mu,lam)) break; nextc+=2; tort++; hare++;} return {det, nextc-det}; }
C Fsim(int mu,int lam){ ll nextc=0; int slow=0,fast=0; while(true){slow++;fast+=2;nextc+=3; if(eqp(slow,fast,mu,lam)) break;} ll det=nextc; slow=0; while(true){ if(eqp(slow,fast,mu,lam)) break; nextc+=2; slow++; fast++;} return {det,nextc-det}; }

int main(){ int M=300,L=300;
  FILE* fp=fopen("/data/build/grid.csv","w");
  fprintf(fp,"mu,lam,An,Bn,Fn,o,Adet,Bdet\n");
  for(int mu=0;mu<=M;mu++) for(int lam=1;lam<=L;lam++){
    RA a=A(mu,lam); C b=Bsim(mu,lam), f=Fsim(mu,lam);
    ll An=a.det+a.ent, Bn=b.n+b.c, Fn=f.n+f.c;
    fprintf(fp,"%d,%d,%lld,%lld,%lld,%d,%lld,%lld\n",mu,lam,An,Bn,Fn,a.o,a.det,b.n);
  }
  fclose(fp);
  printf("done\n");
  return 0; }
