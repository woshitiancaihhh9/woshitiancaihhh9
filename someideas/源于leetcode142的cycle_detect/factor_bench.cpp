#include <bits/stdc++.h>
using namespace std;
using u64=unsigned long long; using u128=__uint128_t; using ll=long long;

// ---------- 64-bit modular arithmetic ----------
static inline u64 mulmod(u64 a,u64 b,u64 m){return (u128)a*b%m;}
static inline u64 powmod(u64 a,u64 e,u64 m){u64 r=1%m;a%=m;while(e){if(e&1)r=mulmod(r,a,m);a=mulmod(a,a,m);e>>=1;}return r;}
bool isprime(u64 n){if(n<2)return false;for(u64 p:{2,3,5,7,11,13,17,19,23,29,31,37}){if(n%p==0)return n==p;}u64 d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}for(u64 a:{2,3,5,7,11,13,17,19,23,29,31,37}){u64 x=powmod(a%n,d,n);if(x==0||x==1||x==n-1)continue;bool ok=false;for(int i=1;i<s;i++){x=mulmod(x,x,n);if(x==n-1){ok=true;break;}}if(!ok)return false;}return true;}
u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}

// ---------- Montgomery (odd modulus) ----------
struct Mont{u64 n,ninv,r2,one; 
  Mont(u64 N){n=N;ninv=N;for(int i=0;i<5;i++)ninv*=2-N*ninv;ninv=0ULL-ninv;r2=(u128)(0ULL-N)%N;r2=mulmod(r2,r2,N);/*approx*/ u128 t=((u128)1<<64)%N; r2=mulmod((u64)t,(u64)t,N); one=(u64)(((u128)1<<64)%N);} 
  inline u64 redc(u128 t)const{u64 m=(u64)t*ninv;u128 u=(t+(u128)m*n)>>64;return u>=n?(u64)(u-n):(u64)u;}
  inline u64 mul(u64 a,u64 b)const{return redc((u128)a*b);}
  inline u64 to(u64 a)const{return mul(a,r2);}
  inline u64 from(u64 a)const{return redc((u128)a);}
  inline u64 add(u64 a,u64 b)const{u64 s=a+b;return s>=n?s-n:s;}
};

// ---------- exact (mu,lam) of x->x^2+c mod p ----------
pair<ll,ll> rho_mulam(u64 p,u64 c,u64 x0){
  auto f=[&](u64 x){return (mulmod(x,x,p)+c)%p;};
  // Brent: find lambda
  u64 power=1,lam=1,tort=x0%p,hare=f(tort);
  while(tort!=hare){ if(power==lam){tort=hare;power<<=1;lam=0;} hare=f(hare); lam++; }
  // find mu
  u64 t=x0%p,h=x0%p; for(u64 i=0;i<lam;i++) h=f(h);
  ll mu=0; while(t!=h){t=f(t);h=f(h);mu++;}
  return {mu,(ll)lam};
}

// ---------- closed forms (validated counts = #successor evals) ----------
int ceil_log2(u64 x){int k=0;while(((u64)1<<k)<x)k++;return k;}
ll Aclosed(ll mu,ll lam){ if(mu==0&&lam==1) return 2; int r0=ceil_log2((u64)mu+1); int rl=ceil_log2((u64)((lam+2)/3)); // ceil(lam/3)
  rl=ceil_log2((u64)((lam+2)/3)); int rmin=max(r0,rl); ll prm=(ll)1<<rmin; ll k=mu+lam; return (lam<=prm)?2*k+1:2*k+prm-1; }
ll Bclosed(ll mu,ll lam){ int kB=ceil_log2((u64)max(lam,mu+1)); return ((ll)1<<kB)-1+lam + (lam+2*mu); }
ll Fclosed(ll mu,ll lam){ ll T= (mu==0)?lam: lam*((mu+lam-1)/lam); return 3*T+2*mu; }

mt19937_64 rng(12345);
u64 randprime(u64 lo,u64 hi){while(true){u64 x=lo+rng()%(hi-lo);x|=1;if(isprime(x))return x;}}

int main(int argc,char**argv){
  vector<int> bits={32,40,48,56,60};
  int K=400;
  printf("=== squaring-count statistics (exact, validated closed forms) ===\n");
  printf("%-6s %8s %8s %8s  %7s %7s  %8s %8s  %8s %8s\n","bits","E[A]","E[B]","E[F]","A/B","A/F","A<=B","A<=F","A<B","A<F");
  // store one representative bitsize instances for timing
  vector<array<ll,3>> timingset; // (A,B,F) for bits=48
  for(int b:bits){ u64 half=(u64)1<<(b/2);
    double sA=0,sB=0,sF=0; ll leB=0,leF=0,ltB=0,ltF=0; ll cnt=0;
    for(int i=0;i<K;i++){ u64 p=randprime(half/2,half), q=randprime(half,half*2); if(p==q)continue; u64 N=p*q; if(p>q)swap(p,q);
      u64 c=1+rng()%(p-1), x0=2+rng()%(p-2);
      auto pr=rho_mulam(p,c,x0); ll mu=pr.first,lam=pr.second; if(lam<=0)continue;
      ll A=Aclosed(mu,lam),B=Bclosed(mu,lam),F=Fclosed(mu,lam);
      sA+=A;sB+=B;sF+=F; leB+=(A<=B);leF+=(A<=F);ltB+=(A<B);ltF+=(A<F);cnt++;
      if(b==48) timingset.push_back({A,B,F});
    }
    printf("%-6d %8.1f %8.1f %8.1f  %7.4f %7.4f  %7.1f%% %7.1f%%  %7.1f%% %7.1f%%\n",
      b,sA/cnt,sB/cnt,sF/cnt, sA/sB, sA/sF, 100.0*leB/cnt,100.0*leF/cnt,100.0*ltB/cnt,100.0*ltF/cnt);
  }

  // ---------- real wall-clock with Montgomery + batched gcd ----------
  // build a real 48-bit semiprime to run the squarings on
  u64 half=(u64)1<<24; u64 P=randprime(half/2,half),Q=randprime(half,half*2); u64 N=P*Q; if((N&1)==0)N++;
  Mont M(N); u64 c=M.to(1+rng()%1000);
  auto run=[&](ll cnt,int extra,int W)->double{ // returns nanoseconds
    u64 x=M.to(2), y=M.to(3), acc=1; int wc=0; double g=0; volatile u64 sink=0;
    auto t0=chrono::high_resolution_clock::now();
    for(ll s=0;s<cnt;s++){ x=M.add(M.mul(x,x),c); for(int e=0;e<extra;e++) x=M.mul(x,x); // expensive-oracle inflation
      u64 d=M.from(x); u64 dy=M.from(y); u64 diff=(d>=dy)?d-dy:d-dy+N; acc=mulmod(acc?acc:1,diff?diff:1,N);
      if(++wc==W){ u64 gg=gcd64(acc,N); sink+=gg; acc=1; wc=0; } }
    auto t1=chrono::high_resolution_clock::now(); sink+=x; (void)sink;
    return chrono::duration<double,nano>(t1-t0).count(); };
  // average per-instance wall-clock over the 48-bit timing set, several extra-cost levels
  printf("\n=== wall-clock (48-bit N, Montgomery, batched gcd W=128), mean +/- std us/instance ===\n");
  printf("%-8s %14s %14s %14s   %8s %8s\n","extra","Floyd","Brent","Proposed","spd/B","spd/F");
  for(int extra: {0,2,8,32,64}){
    vector<double> tA,tB,tF;
    int lim=min((int)timingset.size(),120);
    for(int i=0;i<lim;i++){ ll A=timingset[i][0],B=timingset[i][1],F=timingset[i][2];
      tF.push_back(run(F,extra,128)); tB.push_back(run(B,extra,128)); tA.push_back(run(A,extra,128)); }
    auto ms=[&](vector<double>&v,double&mean,double&sd){mean=0;for(double x:v)mean+=x;mean/=v.size();sd=0;for(double x:v)sd+=(x-mean)*(x-mean);sd=sqrt(sd/v.size());};
    double mA,sdA,mB,sdB,mF,sdF; ms(tA,mA,sdA);ms(tB,mB,sdB);ms(tF,mF,sdF);
    printf("%-8d %7.2f+-%4.2f %7.2f+-%4.2f %7.2f+-%4.2f   %8.3f %8.3f\n",extra,
      mF/1000,sdF/1000,mB/1000,sdB/1000,mA/1000,sdA/1000, mB/mA, mF/mA);
  }
  printf("(us per instance; spd = speedup of Proposed vs Brent / Floyd)\n");
  return 0;
}
