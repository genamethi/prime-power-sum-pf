// Additional modes of scan_q_native. Prefix positions count fractional bits from 1.
#pragma once
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <charconv>
#include <set>
#include <optional>
#include <sstream>
#include <limits>
#include <stdexcept>

namespace prefix_scan {
static U number(const std::string& value) {
    U n{};
    auto result = std::from_chars(value.data(), value.data()+value.size(), n);
    if (result.ec != std::errc{} || result.ptr != value.data()+value.size() || !n)
        throw std::runtime_error("expected a positive integer: " + value);
    return n;
}

static std::mutex output_mutex;

// next() returns an empty optional at end of input. Only job allocation is serialized.
template<class Next, class Work>
static void parallel_jobs(U workers, Next next, Work work) {
    std::mutex dispatch;
    std::exception_ptr failure;
    auto worker = [&] {
        try {
            for (;;) {
                decltype(next()) job;
                {
                    std::lock_guard<std::mutex> lock(dispatch);
                    if (failure || !(job = next())) break;
                }
                work(*job);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(dispatch);
            if (!failure) failure = std::current_exception();
        }
        flint_cleanup();
    };
    std::vector<std::thread> pool;
    try {
        for (U i=1; i<workers; ++i) pool.emplace_back(worker);
    } catch (...) {
        std::lock_guard<std::mutex> lock(dispatch);
        failure = std::current_exception();
    }
    worker();
    for (auto& thread : pool) thread.join();
    if (failure) std::rethrow_exception(failure);
}

static mpz_class prefix(U n, U s, U bits) {
    mpz_class power = 0, root;
    mpz_setbit(power.get_mpz_t(), s+n*bits);
    mpz_root(root.get_mpz_t(), power.get_mpz_t(), n);
    return root;
}

static constexpr auto header = "n2,s,j,m1,m2,M,q,run,required,mod240,q2_divides,exact_powers";
static void candidate(U n,U s,U j,U c,U M,U run,U required,bool accepted,const mpz_class& q) {
    mpz_class modulus=q*q, two=2, residue;
    mpz_powm_ui(residue.get_mpz_t(),two.get_mpz_t(),M,modulus.get_mpz_t());
    bool divides=residue==1, exact=false;
    if (divides) {
        mpz_class difference, small;
        mpz_pow_ui(difference.get_mpz_t(),q.get_mpz_t(),n);
        difference -= mpz_class(1) << (n*j+s);
        difference += mpz_class(1) << c;
        if (difference > 1) {
            U e=mpz_remove(small.get_mpz_t(),difference.get_mpz_t(),q.get_mpz_t());
            exact=small==1 && e>=2 && e%2==0 && e<n/2+1 && n-e>=5 && std::gcd(n-e,6UL)==1;
        }
    }
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cout<<n<<','<<s<<','<<j<<','<<n*j+s<<','<<c<<','<<M<<','<<q<<','<<run
             <<','<<required<<','<<accepted<<','<<divides<<','<<exact<<std::endl;
}

// The least admissible gap gives a necessary run length for every possible e.
static U minimum_gap(U n) {
    for (U e=2*((n-1)/4); e>=2; e-=2)
        if (n-e>=5 && std::gcd(n-e,6UL)==1) return n-e;
    throw std::runtime_error("no admissible smaller exponent");
}

static void band(U n,U s,U lo,U hi,U gap,bool compare) {
    U bits=hi+64;
    mpz_class digits=prefix(n,s,bits);
    // Certify an entire band with one zero in its common interval.
    U overlap=std::min(bits,gap*lo);
    if (!compare && overlap>hi) {
        U zero=mpz_scan0(digits.get_mpz_t(),bits-overlap);
        if (zero<bits-hi) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cerr<<"s="<<s<<" exclusion_zero="<<bits-zero<<'\n';
            return;
        }
    }
    // Read through the first zero after hi so all measured runs are exact.
    while (mpz_scan0(digits.get_mpz_t(),0) >= bits-hi) {
        if (bits>(std::numeric_limits<U>::max()-s)/n/2)
            throw std::runtime_error("prefix too long");
        bits*=2;
        digits=prefix(n,s,bits);
    }
    U run=0;
    while (mpz_tstbit(digits.get_mpz_t(),bits-hi-1-run)) ++run;
    mpz_class initial=digits>>(bits-hi);
    U mod15=mpz_fdiv_ui(initial.get_mpz_t(),15);
    std::array<U,2> count{},best{},at{};
    for (U j=hi;;--j) {
        U offset=bits-j, bit=mpz_tstbit(digits.get_mpz_t(),offset);
        if (!bit && (j%2==s%2)) {
            U c=mpz_scan1(digits.get_mpz_t(),offset)-offset;
            U M=n*j+s-c;
            U q15=(mod15+1)%15;
            bool allowed=c>=4 && M%2==0 && (M%4==0 ? q15==1 : q15==7 || q15==13);
            ++count[allowed];
            if (!at[allowed] || run>best[allowed]) {best[allowed]=run; at[allowed]=j;}
            U required=(gap-1)*j;
            if (run>=required) candidate(n,s,j,c,M,run,required,allowed,(digits>>offset)+1);
        }
        if (j==lo) break;
        mod15=((mod15+15-bit)*8)%15; // inverse of 2 modulo 15
        run=bit ? run+1 : 0;
    }
    std::lock_guard<std::mutex> lock(output_mutex);
    if (compare) for (unsigned b=0;b<2;++b)
        std::cerr<<"s="<<s<<" mod240="<<b<<" positions="<<count[b]<<" longest="<<best[b]
                 <<" j="<<at[b]<<" required_at_j="<<(gap-1)*at[b]<<'\n';
}

static int scan(int argc,char** argv,U workers) {
    if (argc<5 || argc>6 || (argc==6 && std::string(argv[5])!="--compare"))
        throw std::runtime_error("usage: scan_q_native prefix N2 J_MIN J_MAX [--compare]");
    U n=number(argv[2]), lo=number(argv[3]), hi=number(argv[4]);
    U max=std::numeric_limits<U>::max();
    if (n<7 || !(n%2) || lo>hi || hi>max-64 || hi+64>(max-n)/n)
        throw std::runtime_error("need odd N2 >= 7 and a representable inclusive j band");
    U gap=minimum_gap(n), total=n_euler_phi(n), done=0;
    std::cout<<header<<std::endl;
    workers=std::min(workers,total);
    std::cerr<<"[....................] 0/"<<total<<" prefixes, "<<workers<<" workers\n";
    U next_s=1;
    parallel_jobs(workers,[&]() -> std::optional<U> {
        while (next_s<n && std::gcd(next_s,n)!=1) ++next_s;
        if (next_s==n) return {};
        return next_s++;
    },[&](U s) {
        band(n,s,lo,hi,gap,argc==6);
        std::lock_guard<std::mutex> lock(output_mutex);
        ++done;
        U bars=static_cast<U>(20.0L*done/total);
        std::cerr<<"["<<std::string(bars,'#')<<std::string(20-bars,'.')<<"] "
                 <<done<<'/'<<total<<" prefixes complete\n";
    });
    return !std::cout;
}

// FLINT factors both coprime halves, then certifies every prime and the product.
static void factor_part(U M,int sign) {
        mpz_class part=mpz_class(1)<<(sign ? M/2 : M);
        part+=sign ? sign : -1;
        mpz_class product=1;
        fmpz_t z;
        fmpz_init(z);
        fmpz_set_mpz(z,part.get_mpz_t());
        fmpz_factor_t factors;
        fmpz_factor_init(factors);
        { std::lock_guard<std::mutex> lock(output_mutex);
          std::cerr<<"M="<<M<<" factoring "<<(M%2?"2^M-1":sign<0?"2^(M/2)-1":"2^(M/2)+1")<<std::endl; }
        fmpz_factor(factors,z);
        for (slong i=0;i<factors->num;++i) {
            if (!fmpz_is_prime(factors->p+i)) throw std::runtime_error("uncertified factor");
            mpz_class prime,power;
            fmpz_get_mpz(prime.get_mpz_t(),factors->p+i);
            mpz_pow_ui(power.get_mpz_t(),prime.get_mpz_t(),factors->exp[i]);
            product*=power;
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout<<M<<','<<prime<<','<<factors->exp[i]<<std::endl;
        }
        fmpz_factor_clear(factors);
        fmpz_clear(z);
        if (product!=part) throw std::runtime_error("factor product mismatch");
}

template<class Next>
static int factor_jobs(U workers, Next next_M) {
    U pending=0;
    std::set<U> seen, half_done;
    std::cout<<"M,prime,multiplicity"<<std::endl;
    using Part = std::pair<U,int>;
    parallel_jobs(workers,[&]() -> std::optional<Part> {
        if (pending) {
            U M=pending;
            pending=0;
            return Part{M,1};
        }
        U M;
        do { M=next_M(); } while (M && !seen.insert(M).second);
        if (!M) return {};
        if (M%2) return Part{M,0};
        pending=M;
        return Part{M,-1};
    },[&](Part job) {
            auto [M,sign]=job;
            factor_part(M,sign);
            std::lock_guard<std::mutex> lock(output_mutex);
            if (!sign || half_done.erase(M)) std::cerr<<"M="<<M<<" complete\n";
            else half_done.insert(M);
    });
    return !std::cout;
}
static int factor_candidates(U workers) {
    std::string line;
    if (!std::getline(std::cin,line) || line!=header)
        throw std::runtime_error("expected the CSV header from prefix mode on stdin");
    return factor_jobs(workers,[&]() -> U {
        if (!std::getline(std::cin,line)) {
            if (std::cin.bad()) throw std::runtime_error("candidate input error");
            return 0;
        }
        std::istringstream row(line);
        std::string field;
        for (int column=0;column<=5;++column)
            if (!std::getline(row,field,',')) throw std::runtime_error("incomplete candidate row");
        return number(field);
    });
}
} // namespace prefix_scan
