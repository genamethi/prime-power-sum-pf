# Cyclotomic reciprocity calculations

These calculations develop the relative-norm suggestion following [Sharifi, *On cyclotomic polynomials, power residues, and reciprocity laws*](https://www.math.ucla.edu/~sharifi/publ/cycl.pdf). The reciprocity convention is his equation (5). The deductions below are for our equations; they are not applications of his prime-value theorems.

## Local conditions supplied by the equations

Retain

```text
r = ord_q(2),  E = lcm(2,r),  rho = v_q(2^r − 1),
a = (2^E − 1)/q^rho,
Aq^rho = 2^M − 1,  A2^m₂ = q^N − 1,
r | M,  gcd(E,N) = 1,  rho ≥ 2 even.
```

Let s be an odd prime dividing r. Since s divides q−1 and s does not divide N,

```text
v_s(q−1) = v_s(A) = eta_s + v_s(M) ≥ 1 + v_s(r),
eta_s = v_s(2^ord_s(2) − 1) ≥ 1.
```

Here s divides A, so ord_s(2) divides M and lifting applies. No assumption eta_s=1 is needed. Consequently, for m=s^b dividing r,

```text
q ≡ 1 (mod s^(b+1)),
q is an m-th power in Q_s.
```

The second assertion follows from the s-adic logarithm on 1+sZ_s. All m-th norm-residue symbols with q therefore equal 1 at primes above s.

For s=2, q≡1 (mod 8) gives r | (q−1)/2. Thus m=2^b dividing r implies q≡1 (mod 2m). If q≡1 (mod 4m), q is an m-th power in Q_2. Otherwise write q=u^(m/2), where u≡5 (mod 8); adjoining an m-th root of q gives the unramified quadratic extension after adjoining the m-th roots of unity. Its norm-residue character is trivial on units. This includes b=1. The cyclotomic expressions used below are units at 2.

In particular, the remaining equations imply the compact necessary condition

```text
E rad(r) | q−1.
```

## A norm whose factors need not be computed

Choose a prime power m=s^b dividing r. Put

```text
B = 2^(E/m),
H_s = Φ_m(B) = (2^E − 1)/(2^(E/s) − 1),
a_s = H_s/q^rho.
```

The integer H_s is independent of b. We have a_s | a, and B has order m modulo q. In K=Q(ζ_m), let β=B−ζ_m. Then Norm(β)=H_s. There is exactly one prime above q dividing β, with valuation rho: the prime where ζ_m reduces to B. Let J be the ideal obtained from (β) by removing its prime factors above q and s.

At every prime P dividing J, the residue field is F_ell for some ell dividing a_s, and ζ_m reduces to B. Thus ell≡1 (mod m). If

```text
(q/J)_m ≠ 1,
```

some such ell satisfies s | ord_ell(q). This contradicts ord_ell(q) | N and gcd(E,N)=1.

## The calculation at q

All expressions in this paragraph are modulo q. For each i with 1<i<m and gcd(i,m)=1, choose u_i satisfying i*u_i≡1 (mod m). Define

```text
W_m = (-1)^rho · a_s / Φ'_m(B)
      · product_(1<i<m, gcd(i,m)=1) (B − B^i)^u_i,

T_m = W_m^((q−1)/m).
```

Changing the representatives u_i does not change T_m. Under the local conditions above, T_m is the reduction of (q/J)_m at ζ_m=B. Hence

```text
T_m ≠ 1  ⇒  s | ord_ell(q) for some ell | a_s  ⇒  contradiction.
```

**Proof of the formula.** Global reciprocity expresses (q/J)_m as the product of (β,q)_m at the primes above q; the symbols above s are trivial by the preceding local calculation. For m=2 the real symbol is also trivial, since β and q are positive.

At the prime where ζ_m reduces to B^i, with i≠1, β is a unit with residue B−B^i. At the exceptional prime,

```text
β/q^rho ≡ a_s/Φ'_m(B) (mod q).
```

This follows by removing that factor from Φ_m(B)=product_i(B−ζ_m^i). At a prime P above q, the tame symbol (β,q)_m is the ((q−1)/m)-th power of the unit part of (-1)^v_P(β)β. Raising the contribution at the i-th embedding to u_i expresses every contribution using the same root ζ_m. Their product is the stated W_m. No generator for J is required.

Thus the test uses the Mersenne quotient modulo q and cyclotomic expressions modulo q. It does not require the factorization of H_s, or a choice of M or N. The script also computes the same product independently by lifting the roots of unity modulo q^(rho+2).

## Small values of m

For rho even, the first cases simplify to

| m | W_m modulo q | necessary condition |
| --- | --- | --- |
| 2 | a_2 | a_2 is a square modulo q |
| 3 | a_3(2B+1) | a_3(2B+1) is a cube modulo q |
| 4 | −4a_2 | a_2 is a fourth power modulo q |

For m=3, use B²+B+1=0. For m=4, use B²=−1; since q≡1 (mod 8), −4 is a fourth power modulo q.

## Calculations

The [full table](cyclotomic-reciprocity-tables.md) compares the q-adic calculation with direct power-residue symbols at the prime factors. All 23 comparisons satisfying the local condition agree. There are 41 complete-factor calculations and 46 independent comparisons between the derivative formula and the q-adic calculation.

Two checks have rho=2:

* **q=1093, r=E=364, m=2.** Here a_2=(2^182+1)/1093²≡101 (mod 1093), and 101^546≡−1. The complete factorization independently gives the same quadratic symbol.
* **q=3511, r=1755, E=3510, m=3.** Here B=2^1170≡2754, a_3=(B²+B+1)/3511²≡2302, and W_3≡−14. Consequently T_3=(-14)^1170≡756≠1. No factorization of this 705-digit H_3 was used. Independently, 262657 is a primitive Mersenne factor at exponent 27 and 3511^87552≡512≠1 modulo 262657.

These two primes already fail other restrictions from the notes. Their role here is to check the calculation with an even valuation, rather than extrapolate from rho=1.

For **q=274177, r=128, m=128**, the local condition is the boundary case q≡1 (mod 256), q≢1 (mod 512). The calculation gives T_m=8, corresponding to ζ_128³. Directly,

```text
2^64+1 = 274177 · 67280421310721,
274177^((67280421310721−1)/128) ≡ 8 (mod 67280421310721).
```

Products can cancel. For **q=433, m=8**, the individual exponents at ell=17,241,38737 are 6,2,0 modulo 8. Their sum is 0, so T_8=1 despite the two nontrivial symbols. The same q has T_9=27≠1.

## Remaining question

Can all these T_m equal 1 when rho is even and the remaining conditions on q and r hold? The formula makes this a question about the residues of a_s and cyclotomic expressions modulo q, without first factoring a. Proving that some T_m is nontrivial would supply the desired obstruction for all N. No such uniform nonvanishing statement is established here.

For further calculation, simplify W_m using relations between cyclotomic units, starting with m=8,9 and odd primes. A trivial product can be examined through the separate relative norms Norm_(Q(ζ_t)/Q(ζ_m))(2−ζ_t), for m | t | E; the code implements those as well.

```sh
gp -fq -s 64M check_cyclotomic_reciprocity.gp
```

The functions are in [cyclotomic_reciprocity.gp](cyclotomic_reciprocity.gp).
