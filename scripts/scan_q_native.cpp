// Build: g++ -O3 -std=c++17 -pthread scripts/scan_q_native.cpp -lflint -lgmpxx -lgmp -o scan_q_native
// Run: ./scan_q_native MAX_Q THREADS > q-native.csv
#include <flint/ulong_extras.h>
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using U = unsigned long;
using Counts = std::array<U, 6>;
enum class Mode { Both, Candidate, Wieferich };

static U valuation(U x, U p) {
    U v = 0;
    while (x % p == 0) { x /= p; ++v; }
    return v;
}

static U rho_at_least_two(U q) {
    mpz_class modulus = mpz_class(q)*q*q, result, base = 2;
    U rho = 2;
    for (;;) {
        mpz_powm_ui(result.get_mpz_t(), base.get_mpz_t(), q-1, modulus.get_mpz_t());
        if (result != 1) return rho;
        ++rho;
        modulus *= q;
    }
}

static bool repeated_factor(U q) {
    if (q <= 0xffffffffUL) return n_powmod2(2, q-1, q*q) == 1;
    mpz_class modulus = mpz_class(q)*q, result, base = 2;
    mpz_powm_ui(result.get_mpz_t(), base.get_mpz_t(), q-1, modulus.get_mpz_t());
    return result == 1;
}

static void scan(U limit, Mode mode, std::atomic<U>& next, std::mutex& output, Counts& counts) {
    constexpr U chunk = 1UL << 20;
    n_primes_t primes;
    n_primes_init(primes);
    while (true) {
        U lo = next.fetch_add(chunk);
        if (lo > limit) break;
        U hi = std::min(limit, lo + chunk - 1);
        n_primes_jump_after(primes, lo-1);
        for (U q; (q = n_primes_next(primes)) <= hi;) {
            ++counts[0];
            if (q == 2) continue;
            U residue = q % 240;
            bool mod240 = residue == 1 || residue == 97 || residue == 193;
            if (!mod240 && mode != Mode::Wieferich) continue;
            U c = valuation(q-1, 2);
            bool c_ok = c >= 4 && c % 2 == 0;
            if (!c_ok && mode != Mode::Wieferich) continue;
            bool order_ok = (q-1)%6 == 0 && n_powmod2(2, (q-1)/6, q) == 1;
            bool structural = mod240 && c_ok && order_ok;
            if (!structural && mode != Mode::Wieferich) continue;
            if (structural) ++counts[1];
            bool repeated = repeated_factor(q);
            if (!repeated && mode != Mode::Both) continue;
            U rho = repeated ? rho_at_least_two(q) : 1;
            if (repeated) ++counts[2];
            if (rho % 2 == 0) ++counts[3];

            n_factor_t factors;
            n_factor_init(&factors);
            n_factor(&factors, q-1, 1);
            U r = q-1;
            for (int i = 0; i < factors.num; ++i) {
                U ell = factors.p[i];
                while (r % ell == 0 && n_powmod2(2, r/ell, q) == 1) r /= ell;
            }
            U E = std::lcm(2UL, r), rad = 1;
            bool ordinary = true;
            for (int i = 0; i < factors.num; ++i) {
                U ell = factors.p[i];
                if (r % ell == 0) rad *= ell;
                if (ell != 2 && ((q-1)/ell) % r != 0 &&
                    !repeated_factor(ell)) ordinary = false;
            }
            bool rad_ok = (q-1) % (static_cast<__uint128_t>(E)*rad) == 0;
            bool exploratory = structural && ordinary && rad_ok;
            if (exploratory) ++counts[4];
            if (!exploratory && !repeated) continue;

            mpz_class Q;
            mpz_ui_pow_ui(Q.get_mpz_t(), q, rho);
            Q -= 2;
            U height = mpz_sizeinbase(Q.get_mpz_t(), 2);
            std::array<bool, 21> flags = {
                q > 3,q%3 == 1,q%8 == 1,mod240,c_ok,order_ok,repeated,rho % 2 == 0,rad_ok,ordinary,
                E%3 != 0 || q%9 == 1,
                E%3 != 0 || q == 7 || q%7 == 1,
                E%4 != 0 || q == 5 || q%5 == 1,
                E%8 != 0 || q == 17 || q%17 == 1,
                E%12 != 0 || q == 13 || q%13 == 1,
                E%16 != 0 || q == 257 || q%257 == 1,
                E%18 != 0 || q == 19 || q%19 == 1,
                E%4 != 0 || c%4 == 0, q%5 == 1 || (c+rho)%4 == 2,
                E >= 2*height, E >= 2*c
            };
            if (std::all_of(flags.begin(), flags.end(), [](bool v){ return v; })) ++counts[5];
            std::lock_guard<std::mutex> lock(output);
            if (repeated) std::cerr << "q=" << q << " rho=" << rho << '\n';
            std::cout << q;
            for (bool flag : flags) std::cout << ',' << flag;
            std::cout << '\n';
        }
    }
    n_primes_clear(primes);
    flint_cleanup();
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        std::cerr << "usage: scan_q_native MAX_Q THREADS [candidate|wieferich [MIN_Q]]\n";
        return 2;
    }
    U limit = std::stoul(argv[1]), workers = std::stoul(argv[2]);
    if (limit < 2 || limit > (1UL << 60) || workers < 1 || workers > 256) return 2;
    Mode mode = Mode::Both;
    if (argc >= 4) {
        std::string option = argv[3];
        if (option == "candidate") mode = Mode::Candidate;
        else if (option == "wieferich") mode = Mode::Wieferich;
        else return 2;
    }
    U first = argc == 5 ? std::stoul(argv[4]) : 2;
    if (first < 2 || first > limit) return 2;
    std::cout << "q,q_gt_3,mod_3,mod_8,mod_240,c_even_ge_4,r_divides_q_minus_1_over_6,"
                 "rho_ge_2,rho_even,E_rad_r_divides_q_minus_1,ordinary_divisors_q_minus_1,"
                 "E_multiple_3_implies_q_mod_9,E_multiple_3_implies_q_mod_7,"
                 "E_multiple_4_implies_q_mod_5,E_multiple_8_implies_q_mod_17,"
                 "E_multiple_12_implies_q_mod_13,E_multiple_16_implies_q_mod_257,"
                 "E_multiple_18_implies_q_mod_19,E_multiple_4_implies_c_multiple_4,"
                 "mod_5_implies_c_plus_rho,height_Q,height_c\n";
    std::vector<Counts> counts(workers);
    std::vector<std::thread> threads;
    std::mutex output;
    std::atomic<U> next{first};
    for (auto& count : counts)
        threads.emplace_back(scan, limit, mode, std::ref(next), std::ref(output), std::ref(count));
    for (auto& thread : threads) thread.join();
    Counts total{};
    for (auto& count : counts)
        for (std::size_t i = 0; i < total.size(); ++i) total[i] += count[i];
    std::cerr << "primes=" << total[0] << " cheap=" << total[1] << " repeated=" << total[2]
              << " even_rho=" << total[3] << " exploratory=" << total[4]
              << " all_flags=" << total[5] << '\n';
    return !std::cout;
}
