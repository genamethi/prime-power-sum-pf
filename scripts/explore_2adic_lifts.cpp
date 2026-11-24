// Build: make explore_2adic_lifts
// Run: ./build/bin/explore_2adic_lifts > data/2adic-lifts.csv
#include <gmpxx.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

// F(x)=x^e(x^N-1)+2^c. F'(x) is odd at every odd x.
static mpz_class root(unsigned c,unsigned e,unsigned N,unsigned bits) {
    mpz_class x=1;
    for (unsigned known=1;known<bits;known*=2) {
        mpz_class mod=mpz_class(1)<<std::min(bits,2*known);
        mpz_class xe,xn,inverse;
        mpz_powm_ui(xe.get_mpz_t(),x.get_mpz_t(),e-1,mod.get_mpz_t());
        mpz_powm_ui(xn.get_mpz_t(),x.get_mpz_t(),N,mod.get_mpz_t());
        mpz_class derivative=xe*((e+N)*xn-e);
        if (!mpz_invert(inverse.get_mpz_t(),derivative.get_mpz_t(),mod.get_mpz_t())) std::abort();
        mpz_class f=x*xe*(xn-1)+(mpz_class(1)<<c);
        x=(x-f*inverse)%mod;
        if (x<0) x+=mod;
    }
    mpz_class modulus=mpz_class(1)<<bits, xe, xn;
    mpz_powm_ui(xe.get_mpz_t(),x.get_mpz_t(),e,modulus.get_mpz_t());
    mpz_powm_ui(xn.get_mpz_t(),x.get_mpz_t(),N,modulus.get_mpz_t());
    mpz_class check=xe*(xn-1)+(mpz_class(1)<<c);
    if (!mpz_divisible_2exp_p(check.get_mpz_t(),bits)) std::abort();
    return x;
}

int main() {
    constexpr unsigned maxM=8192;
    std::cout<<"c,e,N,windows,top4zero,top8zero,max_leading_zeros,M_at_max,min_excess_bits,M_at_min,size_hits,mod240_hits\n";
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (unsigned c=4;c<=48;c+=2) for (unsigned e=2;e<=48;e+=2) {
        for (unsigned N=std::max(5u,e+1);N*c<maxM;N+=2) {
            if (std::gcd(N,6u)!=1) continue;
            mpz_class x=root(c,e,N,maxM+c);
            unsigned windows=0,top4=0,top8=0,longest=0,bestM=0,atMin=0,hits=0,modHits=0;
            int minExcess=0;
            // q >= 3*2^c+1, hence log2(q)>c+1. Both size inequalities
            // below are deliberately weakened to keep an exponent-only grid.
            for (unsigned M=2;M<=maxM;M+=2) {
                if (N*(c+1)>=M || 2*e*(c+1)>M) continue;
                if (M%4==0 ? (c%4!=0 || N%5==0) : (c+e)%4!=2) continue;
                if (std::gcd(M+c,N+e)!=1) continue;
                unsigned K=M+c;
                // Find the highest set bit below K without copying the residue.
                unsigned len=K;
                while (!mpz_tstbit(x.get_mpz_t(),len-1)) --len;
                unsigned zeros=K-len, cap=M/(2*e)+1;
                ++windows; top4+=zeros>=4; top8+=zeros>=8;
                if (!bestM || zeros>longest) {longest=zeros;bestM=M;}
                if (!atMin || int(len)-int(cap)<int(minExcess)) {
                    minExcess=int(len)-int(cap);atMin=M;
                }
                if (len<=cap) {
                    mpz_class q,power;
                    mpz_fdiv_r_2exp(q.get_mpz_t(),x.get_mpz_t(),K);
                    mpz_pow_ui(power.get_mpz_t(),q.get_mpz_t(),e);
                    if (power>(mpz_class(1)<<(M/2))+1) continue;
                    ++hits;
                    auto residue=mpz_fdiv_ui(q.get_mpz_t(),240);
                    modHits+=M%4==0 ? residue==1 : residue==97 || residue==193;
                }
            }
            if (windows) {
                #pragma omp critical
                std::cout<<c<<','<<e<<','<<N<<','<<windows<<','<<top4<<','<<top8<<','<<longest<<','<<bestM<<','<<int(minExcess)<<','<<atMin<<','<<hits<<','<<modHits<<'\n';
            }
        }
    }
}
