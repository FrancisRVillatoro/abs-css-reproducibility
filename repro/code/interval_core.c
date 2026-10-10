/* Independent directed-binary64 verifier for the absCSS Frobenius recurrence.
 * Build: gcc -std=c11 -O3 -fPIC -shared -frounding-math -ffp-contract=off
 *        -fno-fast-math interval_core.c -lm -o interval_core.so
 * Required: IEC 60559 binary64, gradual underflow, working FE_DOWNWARD.
 * No transcendental operation is performed in the recurrence.
 */
#pragma STDC FENV_ACCESS ON
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef __SSE__
#include <xmmintrin.h>
#endif

typedef struct { double lo,hi; } Iv;
static inline Iv iv(double lo,double hi){Iv a={lo,hi};return a;}
static inline Iv pt(double a){return iv(a,a);}
static inline Iv neg(Iv a){return iv(-a.hi,-a.lo);}
static inline Iv add(Iv a,Iv b){return iv(a.lo+b.lo,-((-a.hi)-b.hi));}
static inline Iv sub(Iv a,Iv b){return add(a,neg(b));}
static inline double mn(double a,double b){return a<b?a:b;}
static inline double mx(double a,double b){return a>b?a:b;}
static inline Iv mul(Iv a,Iv b){
    if(a.lo>=0.){
        if(b.lo>=0.)return iv(a.lo*b.lo,-((-a.hi)*b.hi));
        if(b.hi<=0.)return iv(a.hi*b.lo,-((-a.lo)*b.hi));
        return iv(a.hi*b.lo,-((-a.hi)*b.hi));
    }
    if(a.hi<=0.){
        if(b.lo>=0.)return iv(a.lo*b.hi,-((-a.hi)*b.lo));
        if(b.hi<=0.)return iv(a.hi*b.hi,-((-a.lo)*b.lo));
        return iv(a.lo*b.hi,-((-a.lo)*b.lo));
    }
    if(b.lo>=0.)return iv(a.lo*b.hi,-((-a.hi)*b.hi));
    if(b.hi<=0.)return iv(a.hi*b.lo,-((-a.lo)*b.lo));
    return iv(mn(a.lo*b.hi,a.hi*b.lo),-mn((-a.lo)*b.lo,(-a.hi)*b.hi));
}
static inline Iv invpos(Iv a){return iv(1./a.hi,-((-1.)/a.lo));}
static inline Iv divpos(Iv a,Iv b){return mul(a,invpos(b));}
static inline double abmax(Iv a){return mx(fabs(a.lo),fabs(a.hi));}
static inline Iv absiv(Iv a){return iv(a.lo<=0. && a.hi>=0.?0.:mn(fabs(a.lo),fabs(a.hi)),abmax(a));}
static inline Iv ipow(Iv a,int n){Iv z=pt(1.);while(n){if(n&1)z=mul(z,a);n>>=1;if(n)a=mul(a,a);}return z;}
static inline int enter(void){int old=fegetround();fesetround(FE_DOWNWARD);return old;}
static inline void leave(int old){fesetround(old);}

int iv_environment(void){
    if(FLT_RADIX!=2||DBL_MANT_DIG!=53||DBL_MAX_EXP!=1024)return 0;
#ifdef __SSE__
    if(_mm_getcsr() & ((1u<<15)|(1u<<6)))return 0; /* FTZ or DAZ */
#endif
    int old=enter();volatile double a=1.,b=0x1p-53;double l=a+b;double h=-((-a)-b);
    volatile double c=DBL_MIN,d=0.5;double s=c*d;
    int good=(l==1. && h>1. && s>0. && s<DBL_MIN && fegetround()==FE_DOWNWARD);
    leave(old);return good;
}
void iv_operations(int op,Iv a,Iv b,Iv*out){int old=enter();if(op==0)*out=add(a,b);else if(op==1)*out=sub(a,b);else if(op==2)*out=mul(a,b);else *out=divpos(a,b);leave(old);}

typedef struct {int N;Iv A2,Ah,At,C,D,P,delta;Iv *mat;} Prob;
void* iv_problem(Iv P,int N,const Iv* g,Iv delta){
    if(P.lo<=0.||N<2||!iv_environment())return NULL;
    int old=enter();Prob*b=calloc(1,sizeof(Prob));if(!b){leave(old);return NULL;}
    b->N=N;b->P=P;b->delta=delta;
    Iv pinv=invpos(P),pinv2=mul(pinv,pinv);
    b->A2=add(pt(1.),mul(pt(6.),pinv));b->Ah=mul(pt(8.),pinv2);b->At=sub(pt(1.),b->Ah);b->C=add(b->A2,b->Ah);b->D=mul(pt(4.),pinv2);
    size_t size=(size_t)(N+1)*(N+1);b->mat=calloc(size,sizeof(Iv));
    Iv* dp=calloc(N/2+2,sizeof(Iv));dp[0]=pt(1.);for(int j=1;j<=N/2+1;j++)dp[j]=mul(dp[j-1],delta);
    for(int n=1;n<=N;n++){
        int nm=n-1;Iv dn=add(add(pt((double)nm*nm),mul(pt(6.*nm),pinv)),b->Ah);
        Iv den=mul(pt(n),dn);
        for(int k=1;k<=n;k++){
            Iv t;
            if(k&1){t=neg(mul(b->D,mul(pt(2-k),g[k-1])));t=mul(t,dp[(k-1)/2]);}
            else {int m=n-k;
                t=add(mul(b->A2,mul(pt((double)m*(m-1)),g[k])),mul(b->Ah,mul(pt((double)(1-k)*m),g[k])));
                if(k==2)t=add(t,mul(b->At,pt(m)));
                t=add(t,mul(b->C,g[k-2]));t=neg(mul(t,dp[k/2]));
            }
            b->mat[(size_t)n*(N+1)+k]=divpos(t,den);
        }
    }
    free(dp);leave(old);return b;
}
void iv_free(void*ptr){Prob*b=ptr;if(b){free(b->mat);free(b);}}

/* Radius and tstar enclose exact positive rationals supplied by the caller. */
static Iv major(Prob*b,Iv R,Iv om,int mode){
    Iv pi0=iv(3.1399999999999997,3.1400000000000001);
    Iv q=divpos(mul(mul(R,R),b->delta),mul(pi0,pi0));
    Iv den=sub(pt(1.),q);if(den.lo<=0.)return iv(INFINITY,INFINITY);
    Iv G=divpos(mul(pt(4.),q),den);
    Iv H=add(pt(1.),divpos(mul(mul(pt(4.),q),add(pt(1.),q)),mul(den,den)));
    Iv R2=mul(mul(R,R),b->delta);
    Iv S=mul(absiv(b->A2),G);
    S=add(S,mul(absiv(b->Ah),sub(H,pt(1.))));S=add(S,mul(absiv(b->At),R2));
    S=add(S,mul(absiv(b->C),mul(R2,add(pt(1.),G))));
    S=add(S,mul(mul(absiv(b->D),mul(R,H)),om));return S;
}
static double weight(Iv z,int n,Iv rp){
    double w=(n>=2?(double)n*(n-1):(n==1?1.:1.));
    return mul(pt(abmax(z)),mul(pt(w),rp)).hi;
}
static Iv tailbound(double M,int N,Iv R,Iv tstar){
    Iv rho=divpos(tstar,R);return mul(pt(M),divpos(ipow(rho,N+1),sub(pt(1.),rho)));
}
/* Return q, t q', t^2 q'' as six interval entries (re, im). */
static void evpar(const Iv*a,int N,Iv t,double tail,Iv*out){
    Iv z[6];for(int j=0;j<6;j++)z[j]=pt(0.);
    for(int n=N;n>=0;n--){for(int j=0;j<6;j++)z[j]=mul(z[j],t);
        int k=n&1;z[k]=add(z[k],a[n]);z[2+k]=add(z[2+k],mul(pt(n),a[n]));z[4+k]=add(z[4+k],mul(pt((double)n*(n-1)),a[n]));}
    for(int j=0;j<6;j++)out[j]=add(z[j],iv(-tail,tail));
}
/* mode=0 direct, mode=1 high; om is real frequency or 1 respectively. */
int iv_solve(void*ptr,Iv om,int N,Iv R,Iv tstar,Iv t,Iv*out,double*stats){
    Prob*b=ptr;if(!b||N>b->N)return -1;int old=enter();Iv S=major(b,R,om,0);stats[0]=S.hi;
    if(!(N>S.hi)){leave(old);return 1;}
    Iv*a=calloc(N+1,sizeof(Iv));a[0]=pt(1.);Iv rp=pt(1.);double M=1.;
    for(int n=1;n<=N;n++){
        Iv se=pt(0.),so=pt(0.);Iv*row=b->mat+(size_t)n*(b->N+1);
        for(int k=2;k<=n;k+=2)se=add(se,mul(row[k],a[n-k]));
        for(int k=1;k<=n;k+=2)so=add(so,mul(row[k],a[n-k]));
        so=mul(om,so);a[n]=(n&1)?add(se,so):sub(se,so);
        rp=mul(rp,R);M=mx(M,weight(a[n],n,rp));
    }
    Iv tail=tailbound(M,N,R,tstar);stats[1]=tail.hi;stats[2]=M;
    evpar(a,N,t,tail.hi,out);
    Iv f=divpos(add(mul(out[2],out[0]),mul(out[3],out[1])),t);out[6]=f;
    int good=isfinite(f.lo)&&isfinite(f.hi)&&isfinite(tail.hi)?0:2;
    free(a);leave(old);return good;
}
/* C and E low-frequency system, known X1 coefficient is -sin(t)/2. */
int iv_low(void*ptr,Iv z,int N,Iv R,Iv tstar,Iv t,const Iv*x0,const Iv*x1,Iv*out,double*stats){
    Prob*b=ptr;if(!b||N>b->N)return -1;int old=enter();Iv S=major(b,R,add(pt(1.),z),0);stats[0]=S.hi;if(!(N>S.hi)){leave(old);return 1;}
    Iv*c=calloc(N+1,sizeof(Iv)),*e=calloc(N+1,sizeof(Iv));Iv rp=pt(1.);double M=1.;
    for(int n=1;n<=N;n++){
        Iv sc=pt(0.),se=pt(0.);Iv*row=b->mat+(size_t)n*(b->N+1);
        for(int k=2;k<=n;k+=2){sc=add(sc,mul(row[k],c[n-k]));se=add(se,mul(row[k],e[n-k]));}
        for(int k=1;k<=n;k+=2){sc=sub(sc,mul(row[k],add(x1[n-k],mul(z,e[n-k]))));se=add(se,mul(row[k],c[n-k]));}
        c[n]=sc;e[n]=se;rp=mul(rp,R);
        M=mx(M,mx(mx(weight(c[n],n,rp),weight(e[n],n,rp)),mx(weight(x0[n],n,rp),weight(x1[n],n,rp))));
    }
    double tail=tailbound(M,N,R,tstar).hi;stats[1]=tail;stats[2]=M;
    Iv cs=pt(0.),cp=pt(0.),es=pt(0.),ep=pt(0.);
    for(int n=N;n>=0;n--){cs=add(mul(cs,t),c[n]);cp=add(mul(cp,t),mul(pt(n),c[n]));es=add(mul(es,t),e[n]);ep=add(mul(ep,t),mul(pt(n),e[n]));}
    Iv terr=iv(-tail,tail);cs=add(cs,terr);es=add(es,terr);cp=divpos(add(cp,terr),t);ep=divpos(add(ep,terr),t);
    Iv G=add(sub(add(neg(cs),mul(z,mul(cp,cs))),mul(mul(z,pt(.5)),ep)),mul(mul(z,z),mul(ep,es)));
    out[0]=G;out[1]=cs;out[2]=cp;out[3]=es;out[4]=ep;
    free(c);free(e);leave(old);return isfinite(G.lo)&&isfinite(G.hi)?0:2;
}
/* WKB remainder bound. sqrt_sup[i] must enclose sqrt(1-(i/2048)^6). */
double iv_wkb_bound(Iv P,double rbar,double hbar,const double*sqrt_sup,double*rows){
    int old=enter();Iv p2=mul(P,P),ivp2=invpos(p2);
    Iv AA=add(add(mul(pt(4.),p2),mul(pt(18.),P)),pt(18.));
    Iv EE=add(sub(p2,mul(pt(9.),P)),pt(18.));double bmax=0.;
    for(int j=0;j<2048;j++){
        double wl=j/2048.,wh=(j+1)/2048.;Iv w6=ipow(iv(wl,wh),6);
        double rup=(j==0?rbar:mn(rbar,divpos(pt(hbar),pt(wl)).hi));
        Iv beta1=mul(divpos(pt(2.),pt(9.)),mul(mul(pt(rup),pt(rup)),mul(ivp2,absiv(sub(mul(AA,w6),EE)))));
        Iv beta0=mul(divpos(pt(4.),pt(27.)),mul(ipow(pt(rup),3),mul(ivp2,mul(add(mul(pt(2.),mul(AA,w6)),EE),pt(sqrt_sup[j])))));
        rows[j]=add(beta1,beta0).hi;bmax=mx(bmax,rows[j]);
    }
    leave(old);return bmax;
}
