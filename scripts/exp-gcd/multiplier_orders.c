#define _GNU_SOURCE
// Build: make multiplier_orders
// Run with --help for arguments and output files.
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef __linux__
#include <sched.h>
#endif

enum { MAX_N = 2519, MAX_H = 630, MAX_FACTORS = 128, MAX_SLOTS = 4,
       MAX_WIDTH = 32, MAX_ROWS = 850 };
static int *spf, *order2;

static void fail(const char *msg) { fprintf(stderr, "%s\n", msg); exit(EXIT_FAILURE); }
static void *allocate(size_t count, size_t size) {
    void *p = calloc(count, size);
    if (!p) fail("out of memory");
    return p;
}
static int gcd(int a, int b) {
    while (b) { int rem = a % b; a = b; b = rem; }
    return a;
}
static int power_mod(int a, int e, int m) {
    int64_t x = a, result = 1;
    while (e) {
        if (e % 2) result = result * x % m;
        x = x * x % m;
        e /= 2;
    }
    return (int)result;
}
static void init_orders(int limit) {
    spf = allocate(limit + 1, sizeof *spf);
    order2 = allocate(limit + 1, sizeof *order2);
    for (int p = 2; p <= limit; ++p) {
        if (spf[p]) continue;
        spf[p] = p;
        if (p <= limit/p)
            for (int j = p*p; j <= limit; j += p)
                if (!spf[j]) spf[j] = p;
        if (p == 2) continue;
        int order = p-1, rest = p-1;
        while (rest > 1) {
            int f = spf[rest];
            while (order % f == 0 && power_mod(2, order/f, p) == 1) order /= f;
            do { rest /= f; } while (rest % f == 0);
        }
        order2[p] = order;
    }
}

typedef struct { int p, power; } Factor;
typedef struct { int count; Factor list[MAX_FACTORS]; } Factors;

// Store an lcm by its maximal prime powers. Its integer value need never fit a word.
// Descending order permits a single pass when ord_p(2) introduces smaller primes.
static void require_power(Factors *f, int p, int power) {
    int i = 0;
    while (i < f->count && f->list[i].p > p) ++i;
    if (i < f->count && f->list[i].p == p) {
        if (f->list[i].power < power) f->list[i].power = power;
        return;
    }
    if (f->count == MAX_FACTORS) fail("increase MAX_FACTORS");
    memmove(f->list+i+1, f->list+i, (size_t)(f->count-i)*sizeof(Factor));
    f->list[i] = (Factor){p,power};
    ++f->count;
}
static void require_divisor(Factors *f, int x) {
    while (x > 1) {
        int p = spf[x], power = 1;
        do { x /= p; power *= p; } while (x % p == 0);
        require_power(f,p,power);
    }
}
static int divides_required(const Factors *f, int divisor) {
    for (int i = 0; i < f->count && divisor > 1; ++i)
        divisor /= gcd(divisor, f->list[i].power);
    return divisor == 1;
}
static int contains(const int *a, int n, int x) {
    for (int i = 0; i < n; ++i) if (a[i] == x) return 1;
    return 0;
}

typedef struct {
    int t, w;                 // w = t^j, exact order of q mod the witness
    int tested, passed, kept, prime[MAX_WIDTH];
    int *all;
} Slot;
typedef struct {
    int n, shared, hmax, forbidden[MAX_SLOTS], nforbidden;
    int hcount, allowed[MAX_H+1], hlist[MAX_H], slots;
    Slot slot[MAX_SLOTS];
    int *minimal, nminimal, excluded;
    uint64_t tuples, repeated, passed;
    int good[MAX_SLOTS], good_h, bad[MAX_SLOTS], have_bad;
} Gap;

static int forbidden_divides(const Gap *g, const Factors *required) {
    for (int i=0; i<g->nforbidden; ++i)
        if (divides_required(required,g->forbidden[i])) return 1;
    return 0;
}

// For one h, close the necessary divisibilities of M under P|E => ord_P(2)|M.
static int compatible(const Gap *g, int h, const Factors *input, const int *witness, int nw) {
    Factors required = *input;
    require_divisor(&required, 2*h);            // 2|E
    for (int i=0; i<required.count; ++i) {
        Factor f = required.list[i];
        if (f.p == 2 || h % f.power == 0) continue;
        // The full power in M cannot fit in h, so P must divide E.
        if (g->n % f.p == 0 || contains(witness,nw,f.p)) return 0;
        require_divisor(&required,order2[f.p]);
    }
    return !forbidden_divides(g,&required);
}
static Factors base_required(const Gap *g) {
    Factors f = {0};
    require_divisor(&f, 2);
    if (g->shared) require_divisor(&f,order2[1093]);
    return f;
}

// Any P ineligible for E has to occur with full multiplicity in h.
// This inexpensive lower divisor skips impossible h before the recursive order test.
static int first_multiplier(const Gap *g, const Factors *f, const int *witness, int nw) {
    if (forbidden_divides(g,f)) return 0;
    int needed = 1;
    for (int i = 0; i < f->count; ++i) {
        Factor x = f->list[i];
        if (x.p == 2) continue;
        int outside = g->n % x.p == 0 || contains(witness,nw,x.p);
        if (!outside) {
            int order = order2[x.p];
            int even_order = order % 2 ? 2*order : order;
            for (int j=0; j<g->nforbidden; ++j)
                if (even_order % g->forbidden[j] == 0) { outside = 1; break; }
        }
        if (outside) {
            if (needed > g->hmax/x.power) return 0;
            needed *= x.power;
        }
    }
    if (needed == 1) {
        for (int i = 0; i < g->hcount; ++i)
            if (compatible(g,g->hlist[i],f,witness,nw)) return g->hlist[i];
    } else {
        for (int h = needed; h <= g->hmax; h += needed)
            if (g->allowed[h] && compatible(g,h,f,witness,nw)) return h;
    }
    return 0;
}
static void prepare_gap(Gap *g) {
    int n = g->n, rest = n;
    g->hmax = (n+1)/4;                          // strict h < (N+2)/4
    while (rest > 1) {
        int t = spf[rest], w = 1;
        if (!(t == 1093 && g->shared) && gcd(n,t-1)==1)
            g->forbidden[g->nforbidden++] = order2[t];
        do {
            w *= t;
            if (g->slots == MAX_SLOTS) fail("increase MAX_SLOTS");
            g->slot[g->slots++] = (Slot){.t=t, .w=w};
            rest /= t;
        } while (rest % t == 0);
    }
    Factors f = base_required(g);
    for (int h=1; h<=g->hmax; ++h) {
        int common = gcd(n,h);                   // gcd(E,N)=1 gives gcd(h,N)=gcd(M,N)
        if (gcd(common,3570)!=1 || (common>1 && n<5*common)) continue;
        if (!compatible(g,h,&f,NULL,0)) continue;
        g->allowed[h] = 1;
        g->hlist[g->hcount++] = h;
    }
}
static void scan_divisors(Gap *g, int bound) {
    unsigned char *marked = allocate(bound+1, sizeof *marked);
    g->minimal = allocate(bound/2+1,sizeof *g->minimal);
    for (int b=1; b<=bound; b+=2) {
        if (!marked[b]) {
            Factors f = base_required(g);
            require_divisor(&f,b);
            if (first_multiplier(g,&f,NULL,0)) continue;
            g->minimal[g->nminimal++] = b;
            for (int j=b; j<=bound; j+=2*b) marked[j]=1;
        }
        ++g->excluded;
    }
    free(marked);
}
static void scan_slot(Gap *g, Slot *s, int limit, int width, int keep_all) {
    // Phi_(t^j)(q) = t mod 16. If t divides q-1, remove its single t first.
    // The remaining factor still exceeds 1 and has only primes of exact order t^j.
    int modulus = 1;
    if (keep_all) s->all=allocate(limit/(2*s->w)+1,sizeof *s->all);
    if (!g->shared && s->t % 16 != 1)
        modulus = 2*((s->t-1) & -(s->t-1));
    for (int p=2*s->w+1; p<=limit; p+=2*s->w) {
        if (spf[p] != p || (modulus>1 && p%modulus==1)) continue;
        ++s->tested;
        Factors f = base_required(g);
        require_divisor(&f,order2[p]);
        if (!first_multiplier(g,&f,&p,1)) continue;
        if (s->all) s->all[s->passed]=p;
        ++s->passed;
        if (s->kept < width) s->prime[s->kept++] = p;
    }
}
static void test_tuple(Gap *g, int pos, int *primes) {
    if (pos == g->slots) {
        ++g->tuples;
        Factors f = base_required(g);
        for (int i=0; i<pos; ++i) {
            if (contains(primes,i,primes[i])) { ++g->repeated; return; }
            require_divisor(&f,order2[primes[i]]);
        }
        int h = first_multiplier(g,&f,primes,pos);
        if (h) {
            ++g->passed;
            if (!g->good_h) {
                memcpy(g->good,primes,(size_t)pos*sizeof *primes);
                g->good_h=h;
            }
        } else if (!g->have_bad) {
            memcpy(g->bad,primes,(size_t)pos*sizeof *primes);
            g->have_bad=1;
        }
        return;
    }
    if (pos >= MAX_SLOTS) fail("invalid tuple depth");
    Slot *s = &g->slot[pos];
    for (int i=0; i<s->kept; ++i) {
        primes[pos]=s->prime[i];
        test_tuple(g,pos+1,primes);
    }
}

static FILE *output(const char *dir, const char *name) {
    char path[4096];
    if (snprintf(path,sizeof path,"%s/%s",dir,name) >= (int)sizeof path) fail("output path too long");
    FILE *f = fopen(path,"w");
    if (!f) { perror(path); exit(EXIT_FAILURE); }
    return f;
}
static void finish(FILE *f) { if (ferror(f) || fclose(f)) fail("writing results failed"); }
static void numbers(FILE *f, const int *a, int n) {
    for (int i=0; i<n; ++i) fprintf(f,"%s%d",i?" ":"",a[i]);
}
static void ranges(FILE *f, const int *a, int n) {
    for (int i=0; i<n; ++i) {
        int start=i;
        while (i+1<n && a[i+1]==a[i]+1) ++i;
        fprintf(f,"%s%d",start?" ":"",a[start]);
        if (i!=start) fprintf(f,"-%d",a[i]);
    }
}
static void gaps_csv(const char *dir, Gap *rows, int count) {
    FILE *g=output(dir,"gaps.csv"), *d=output(dir,"forbidden-divisors.csv"), *w=output(dir,"orders.csv");
    fputs("N,shared_with_q_minus_1,forbidden_orders,h_max,h_count,h_ranges,odd_divisors_excluded,minimal_divisor_count,tuples_checked,repeated_prime_tuples,distinct_tuples_excluded,tuples_passing,first_passing_tuple,first_passing_h,first_excluded_distinct_tuple\n",g);
    fputs("N,shared_with_q_minus_1,minimal_forbidden_odd_divisors\n",d);
    fputs("N,shared_with_q_minus_1,required_q_order,local_primes,single_primes_passing,primes_used_in_tuples\n",w);
    for (int i=0; i<count; ++i) {
        Gap *r = &rows[i]; int shared = r->shared ? 1093 : 1;
        fprintf(g,"%d,%d,",r->n,shared); numbers(g,r->forbidden,r->nforbidden);
        fprintf(g,",%d,%d,",r->hmax,r->hcount); ranges(g,r->hlist,r->hcount);
        fprintf(g,",%d,%d,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",",
                r->excluded,r->nminimal,r->tuples,r->repeated,r->tuples-r->repeated-r->passed,r->passed);
        if (r->good_h) numbers(g,r->good,r->slots);
        fputc(',',g); if (r->good_h) fprintf(g,"%d",r->good_h);
        fputc(',',g); if (r->have_bad) numbers(g,r->bad,r->slots); fputc('\n',g);
        fprintf(d,"%d,%d,",r->n,shared); numbers(d,r->minimal,r->nminimal); fputc('\n',d);
        for (int j=0; j<r->slots; ++j) {
            Slot *s=&r->slot[j];
            fprintf(w,"%d,%d,%d,%d,%d,",r->n,shared,s->w,s->tested,s->passed);
            numbers(w,s->prime,s->kept); fputc('\n',w);
        }
    }
    finish(g); finish(d); finish(w);
}
typedef struct { int b, count; } Frequency;
typedef struct { int count, partner, h; } PairRow;
typedef struct {
    Gap *gap;
    PairRow *data;
    atomic_int next;
} PairWork;

// Exhaust every combination of the two individually viable prime lists for one N.
static void *pair_worker(void *arg) {
    PairWork *job=arg;
    Gap *g=job->gap;
    Slot *a=&g->slot[0], *b=&g->slot[1];
    for (;;) {
        int i=atomic_fetch_add_explicit(&job->next,1,memory_order_relaxed);
        if (i>=a->passed) break;
        Factors first=base_required(g);
        require_divisor(&first,order2[a->all[i]]);
        for (int j=0; j<b->passed; ++j) {
            int pair[]={a->all[i],b->all[j]};
            if (pair[0]==pair[1]) continue;
            Factors required=first;
            require_divisor(&required,order2[pair[1]]);
            int h=first_multiplier(g,&required,pair,2);
            if (!h) continue;
            PairRow *r=&job->data[i];
            if (!r->count) { r->partner=pair[1]; r->h=h; }
            ++r->count;
        }
    }
    return NULL;
}
static void run_workers(int threads, void *(*worker)(void *), void *work) {
    pthread_t *workers=allocate(threads,sizeof *workers);
    for (int i=0; i<threads; ++i) {
        int err=pthread_create(&workers[i],NULL,worker,work);
        if (err) fail(strerror(err));
    }
    for (int i=0; i<threads; ++i) {
        int err=pthread_join(workers[i],NULL);
        if (err) fail(strerror(err));
    }
    free(workers);
}
static uint64_t scan_all_pairs(const char *dir, Gap *g, int threads) {
    PairRow *data=allocate(g->slot[0].passed+1,sizeof *data);
    PairWork job={.gap=g,.data=data,.next=0};
    run_workers(threads,pair_worker,&job);
    FILE *f=output(dir,"all-pairs.csv");
    fputs("N,first_order,first_prime,second_order,individually_passing_second_primes,passing_partners,first_passing_partner,first_passing_h\n",f);
    uint64_t passed=0;
    for (int i=0; i<g->slot[0].passed; ++i) {
        PairRow r=data[i]; passed+=r.count;
        fprintf(f,"%d,%d,%d,%d,%d,%d,",g->n,g->slot[0].w,g->slot[0].all[i],
                g->slot[1].w,g->slot[1].passed,r.count);
        if (r.count) fprintf(f,"%d,%d",r.partner,r.h); else fputc(',',f);
        fputc('\n',f);
    }
    finish(f);
    free(data);
    return passed;
}
static int most_frequent(const void *x, const void *y) {
    const Frequency *a=x,*b=y;
    return a->count==b->count ? a->b-b->b : b->count-a->count;
}
static void index_md(const char *dir, Gap *rows, int count, int bmax, int pmax, int width, int threads,
                     const Gap *focus, uint64_t deep_pass) {
    FILE *f=output(dir,"README.md");
    int htotal=0, hpass=0, no_h=0, empty_pool=0, no_joint=0, cases_passing=0;
    uint64_t families=0, tuples=0, passed=0, repeated=0;
    Frequency *freq=allocate(bmax/2+1,sizeof *freq);
    for (int j=0; j<=(bmax-1)/2; ++j) freq[j].b=2*j+1;
    for (int i=0; i<count; ++i) {
        Gap *g=&rows[i];
        htotal+=g->hmax; hpass+=g->hcount; no_h+=!g->hcount;
        families+=g->excluded; tuples+=g->tuples; passed+=g->passed; repeated+=g->repeated;
        if (!g->tuples) ++empty_pool;
        else if (!g->passed) ++no_joint;
        else ++cases_passing;
        for (int j=0; j<g->nminimal; ++j) ++freq[g->minimal[j]/2].count;
    }
    qsort(freq,(bmax+1)/2,sizeof *freq,most_frequent);
    fputs("# Multiplier and order restrictions\n\n",f);
    fprintf(f,"All %d odd gaps in `5 ≤ N ≤ 2519` with `3 ∤ N` are included. "
              "The two possibilities for `1093 | q−1` at `N = 1093` give %d rows. "
              "The hypotheses are the remaining case of the notes, with `M = hE`, "
              "`E | q−1` and `gcd(E,N) = 1`.\n\n",count-1,count);
    fputs("| Results | Contents |\n|---|---|\n"
          "| [Gaps and multipliers](gaps.csv) | One row per gap and branch, permitted multiplier ranges, assignment counts and examples |\n"
          "| [Forbidden divisors](forbidden-divisors.csv) | Minimal odd divisors proved impossible in M, indexed by N |\n"
          "| [Candidate primes by order](orders.csv) | Counts through the prime limit and the actual primes used in the joint tests |\n",f);
    fprintf(f,"| [All pairs for N = %d](all-pairs.csv) | Count and first example of jointly passing partners for each individually passing prime of order %d |\n\n",
            focus->n,focus->slot[0].w);
    fputs("## Run totals\n\n| Calculation | Result |\n|---|---:|\n",f);
    fprintf(f,"| Multipliers excluded out of those in the height range | %d / %d |\n"
              "| Branches with every multiplier excluded | %d |\n"
              "| Proved exclusions of an odd divisor of M, including multiples, out of divisor and branch combinations | %"PRIu64" / %"PRIu64" |\n"
              "| Joint assignments tested | %"PRIu64" |\n"
              "| Assignments reusing a prime for distinct orders | %"PRIu64" |\n"
              "| Distinct assignments excluded by the multiplier tests | %"PRIu64" |\n"
              "| Assignments passing those tests | %"PRIu64" |\n"
              "| Branches missing a candidate list within the prime limit | %d |\n"
              "| Branches with lists but no passing joint assignment within the tuple cutoff | %d |\n"
              "| Branches with a passing joint assignment | %d |\n\n",
              htotal-hpass,htotal,no_h,families,(uint64_t)count*((bmax+1)/2),tuples,repeated,
              tuples-repeated-passed,passed,empty_pool,no_joint,cases_passing);
    fputs("## Most frequent minimal forbidden divisors\n\n"
          "| Odd divisor | Gap and branch rows forbidding it minimally |\n|---:|---:|\n",f);
    for (int i=0; i<(bmax+1)/2 && i<15 && freq[i].count; ++i)
        fprintf(f,"| %d | %d |\n",freq[i].b,freq[i].count);
    free(freq);
    fprintf(f,"\n## Reading the results\n\n"
              "- **Divisibility.** For each odd integer B through %d the program assumes B divides M and tries every permitted h. If all fail, every multiple of B is excluded, even beyond the cutoff. Only minimal B under divisibility are printed. Within the tested range, B is excluded exactly when a listed divisor divides B. Testing divisors allows any odd part of M that is a multiple of B. The earlier exclusion of odd parts exactly 1 and 3 is a separate result.\n"
              "- **Primes.** For each prime power w = t^j dividing N, the program tries every prime ℓ ≤ %d with w dividing ℓ−1. Since Φ_w(q) divides A, it needs a prime of exact q-order w. If t does not divide q−1, then Φ_w(q) ≡ t modulo 16. For t ≢ 1 modulo 16 at least one such prime differs from 1 modulo the smallest power of two distinguishing t from 1. This uses only c ≥ 4. In the shared 1093 branch, dividing Φ_1093(q) by 1093 leaves a factor greater than 1 with primes of order 1093.\n"
              "- **Joint assignments.** For each w take the first %d individually passing primes in increasing order, or all if fewer. Test their full Cartesian product. Primes must be distinct and all orders must work with one h. Tuple entries follow the order of the rows in [orders.csv](orders.csv). A passing tuple is a set of necessary local conditions. Failure of one tuple proves that assignment impossible. Exhausting bounded lists cannot exclude N.\n"
              "- **Multiplier ranges.** h starts at 1 and satisfies h < (N+2)/4. With g = gcd(h,N), require gcd(g,3570) = 1 and N ≥ 5g for g > 1. For t | N coprime to q(q−1) with gcd(t−1,N) = 1, ord_t(2) cannot divide M. The order table verifies that 1093 is the only odd base-2 Wieferich prime through 2519 and that it cannot equal q because q ≡ 1 modulo 48. This accounts for all branches of gcd(N,q−1).\n"
              "- **Divisibility propagation.** Start with 2h dividing M, any proposed B and the binary orders of the selected primes dividing M. In the shared branch also require ord_1093(2) dividing M. If a required power P^k for odd P does not divide h, then P divides E and q−1. Discard the assignment if P divides N or was assigned nontrivial q-order. Otherwise require ord_P(2) dividing M and repeat. Its prime factors are smaller than P. Discard if any forbidden order now divides M. Every row uses these rules; further conditions from the notes or from a particular q can give additional exclusions.\n",bmax,pmax,width);
    for (int i=0; i<count; ++i) if (rows[i].n==1793) {
        Gap *r=&rows[i];
        fprintf(f,"\n## N = 1793\n\nThe forbidden orders are 10 and 162, so 5 and 81 cannot divide the even gap M. "
                  "The program excludes %d of the %d odd divisors tested. ",r->excluded,(bmax+1)/2);
        if (bmax>=451) fputs("The divisor 451 = 11·41 is excluded by the common multiplier bound. ",f);
        fprintf(f,"Among the tested tuples, %"PRIu64" distinct assignments fail and %"PRIu64" pass.\n",
                r->tuples-r->repeated-r->passed,r->passed);
    }
    fprintf(f,"\nFor **N = %d**, the [all-pairs check](all-pairs.csv) also exhausted the individually passing "
              "primes of orders %d and %d through %d, without the per-order tuple cutoff. "
              "Of %d × %d = %"PRIu64" combinations, **%"PRIu64"** admit a common multiplier under the stated tests.\n",
            focus->n,focus->slot[0].w,focus->slot[1].w,pmax,focus->slot[0].passed,focus->slot[1].passed,
            (uint64_t)focus->slot[0].passed*focus->slot[1].passed,deep_pass);
    fprintf(f,"\n## Reproduce\n\nBuild the [C program](../../scripts/exp-gcd/multiplier_orders.c) and repeat this run from the repository root\n\n"
              "```sh\nmake multiplier_orders\n"
              "./build/bin/multiplier_orders data/multiplier-orders %d %d %d %d %d\n```\n\n"
              "Arguments after the directory are the odd divisor limit, prime limit and number of primes per order in joint tests. "
              "The next argument sets the POSIX worker count. Its default is the number of available CPUs (process affinity on Linux). "
              "The last argument selects the gap for the all-pairs check; it must have exactly two prime factors counted with multiplicity. "
              "Rerunning replaces the five named result files.\n",bmax,pmax,width,threads,focus->n);
    finish(f);
}
static int arg_number(const char *s, int lo, int hi) {
    char *end;
    errno=0;
    long v=strtol(s,&end,10);
    if (errno || end==s || *end || v<lo || v>hi) fail("argument outside supported range; see --help");
    return (int)v;
}
typedef struct {
    Gap *rows;
    int count, bmax, pmax, width, focus;
    atomic_int next;
} Work;
static void *scan_worker(void *arg) {
    Work *w = arg;
    for (;;) {
        int i=atomic_fetch_add_explicit(&w->next,1,memory_order_relaxed);
        if (i>=w->count) break;
        Gap *g=&w->rows[i];
        prepare_gap(g);
        scan_divisors(g,w->bmax);
        for (int j=0; j<g->slots; ++j) scan_slot(g,&g->slot[j],w->pmax,w->width,g->n==w->focus);
        int primes[MAX_SLOTS];
        test_tuple(g,0,primes);
    }
    return NULL;
}
static int available_cpus(void) {
#ifdef __linux__
    cpu_set_t set;
    if (sched_getaffinity(0,sizeof set,&set)==0 && CPU_COUNT(&set)>0)
        return CPU_COUNT(&set);
#endif
    long n=sysconf(_SC_NPROCESSORS_ONLN);
    return n>0 && n<65536 ? (int)n : 1;
}
int main(int argc, char **argv) {
    if (argc<2 || argc>7 || !strcmp(argv[1],"--help")) {
        puts("Usage: multiplier_orders OUTPUT_DIR [ODD_DIVISOR_LIMIT [PRIME_LIMIT [TUPLE_WIDTH [THREADS [PAIR_GAP]]]]]\n"
             "Defaults: 4095, 1000000, 12. Limits: 1..100001, 3..10000000, 1..32.\n"
             "THREADS defaults to available CPUs (process affinity on Linux).\n"
             "PAIR_GAP defaults to 1793; it must have two prime factors counted with multiplicity.\n"
             "Tests all admissible N through 2519. Writes README.md, gaps.csv, orders.csv, all-pairs.csv and\n"
             "forbidden-divisors.csv. The parent of OUTPUT_DIR must exist.");
        return (argc>=2 && !strcmp(argv[1],"--help")) ? 0 : 2;
    }
    int bmax=argc>2?arg_number(argv[2],1,100001):4095;
    int pmax=argc>3?arg_number(argv[3],3,10000000):1000000;
    int width=argc>4?arg_number(argv[4],1,MAX_WIDTH):12;
    int threads=argc>5?arg_number(argv[5],1,65535):available_cpus();
    int selected=argc>6?arg_number(argv[6],5,MAX_N):1793;
    if (bmax%2==0) --bmax;
    if (mkdir(argv[1],0777) && errno!=EEXIST) { perror(argv[1]); return 1; }
    int limit=bmax>pmax?bmax:pmax;
    init_orders(limit>MAX_N?limit:MAX_N);
    for (int p=3; p<=MAX_N; p+=2) if (spf[p]==p) {
        int exceptional = power_mod(2,p-1,p*p)==1;
        if (exceptional != (p==1093) || (exceptional && p%48==1))
            fail("recheck the branches for q dividing N or sharing a factor with N");
    }
    Gap *rows=allocate(MAX_ROWS,sizeof *rows);
    int count=0;
    for (int n=5; n<=MAX_N; n+=2) {
        if (n%3==0) continue;
        for (int branch=0; branch<= (n==1093); ++branch) {
            if (count==MAX_ROWS) fail("increase MAX_ROWS");
            rows[count].n=n;
            rows[count++].shared=branch;
        }
    }
    if (threads>count) threads=count;
    Gap *focus=NULL;
    for (int i=0; i<count; ++i) if (rows[i].n==selected) focus=&rows[i];
    int first=spf[selected], other=selected/first;
    if (!focus || other==1 || spf[other]!=other) fail("PAIR_GAP must be an admissible product of two primes");
    Work work={.rows=rows,.count=count,.bmax=bmax,.pmax=pmax,.width=width,.focus=selected,.next=0};
    fprintf(stderr,"Starting %d POSIX workers for %d gap/branch rows\n",threads,count);
    run_workers(threads,scan_worker,&work);
    fprintf(stderr,"Checking all %d x %d individually passing prime pairs for N=%d\n",
            focus->slot[0].passed,focus->slot[1].passed,selected);
    uint64_t deep_pass=scan_all_pairs(argv[1],focus,threads);
    gaps_csv(argv[1],rows,count);
    index_md(argv[1],rows,count,bmax,pmax,width,threads,focus,deep_pass);
    free(focus->slot[0].all); free(focus->slot[1].all);
    for (int i=0; i<count; ++i) free(rows[i].minimal);
    free(rows); free(spf); free(order2);
    printf("Wrote %d gap/branch rows to %s\n",count,argv[1]);
    return 0;
}
