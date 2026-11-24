# Comparing q in Sage

Start `pixi run sage` from the repository root and load both scripts.

```python
load("scripts/inspect_reciprocity.sage")
load("scripts/q_comparison_report.sage")
from itertools import islice
```

## Check conditions and select primes

Pass one prime, a list of primes, or a finite iterator to `q_conditions`.
It returns a list of Boolean rows in input order, even for a single
prime. Pair each row with `Q_CONDITIONS` to see what the flags mean.

`q_conditions` evaluates all the columns in the chosen mode. Choose
which conditions matter when you filter its output. In the example
below, `True` requires a pass, `False` requires a failure, and omitting
a name leaves that condition unrestricted.

```python
qs = list(islice(Primes(), 1000))
rows = q_conditions(qs, silent=True, workers=4, csv_path="q-conditions.csv")

required = {
    "mod_240": True,                    # q ≡ 1, 97, or 193 mod 240
    "r_divides_q_minus_1_over_6": True,  # 6r divides q−1
    "rho_ge_2": False,                  # ρ = 1 for these odd primes
}
flags_by_q = {
    q: dict(zip(Q_CONDITIONS, row))
    for q, row in zip(qs, rows)
}
keep = [
    q for q, flags in flags_by_q.items()
    if all(flags[name] == expected for name, expected in required.items())
]
```

To inspect one row, use `flags_by_q[433]`. For a prime outside
the batch, use `dict(zip(Q_CONDITIONS, q_conditions(q, silent=True)[0]))`.
The example deliberately keeps ρ = 1 so there are more examples for
reciprocity. To require every checked condition, use
`keep = [q for q, row in zip(qs, rows) if all(row)]`.

Omit `silent=True` to print the table. Set `csv_path` to save it
regardless of printing. With `cheap=True`, use `Q_CHEAP_CONDITIONS`
instead of `Q_CONDITIONS` when pairing names with rows. The cheap
checks use congruences and modular powers. The full check factors
`q−1` to determine r. `workers` sets the number of processes and
`batch_size` sets how many primes go in each batch.

### Condition reference

Here r is the order of 2 modulo q, E = lcm(2,r), c = v₂(q−1),
ρ = v_q(2ʳ−1), and Q = q^ρ. An implication returns `True`
when its hypothesis does not hold. The forced congruences modulo a
prime ℓ apply only when q ≠ ℓ. For q = 2 every flag is `False`, since
r is undefined.

| Column | What it checks | Mode |
| --- | --- | --- |
| `q_gt_3` | q > 3 | Cheap or full |
| `mod_3` | q ≡ 1 mod 3 | Cheap or full |
| `mod_8` | q ≡ 1 mod 8 | Cheap or full |
| `mod_240` | q ≡ 1, 97, or 193 mod 240 | Cheap or full |
| `c_even_ge_4` | c is even and at least 4 | Cheap or full |
| `r_divides_q_minus_1_over_6` | 6r divides q−1, checked by modular exponentiation | Cheap or full |
| `rho_ge_2` | ρ ≥ 2 | Cheap or full |
| `rho_even` | ρ is even | Cheap or full |
| `E_rad_r_divides_q_minus_1` | E · rad(r) divides q−1, where rad(r) is the product of the distinct primes dividing r | Full |
| `ordinary_divisors_q_minus_1` | For each odd prime ℓ dividing q−1 with ℓ² ∤ 2^(ℓ−1)−1, r divides (q−1)/ℓ | Full |
| `E_multiple_3_implies_q_mod_9` | If 3 divides E, then q ≡ 1 mod 9 | Full |
| `E_multiple_3_implies_q_mod_7` | If 3 divides E and q ≠ 7, then q ≡ 1 mod 7 | Full |
| `E_multiple_4_implies_q_mod_5` | If 4 divides E and q ≠ 5, then q ≡ 1 mod 5 | Full |
| `E_multiple_8_implies_q_mod_17` | If 8 divides E and q ≠ 17, then q ≡ 1 mod 17 | Full |
| `E_multiple_12_implies_q_mod_13` | If 12 divides E and q ≠ 13, then q ≡ 1 mod 13 | Full |
| `E_multiple_16_implies_q_mod_257` | If 16 divides E and q ≠ 257, then q ≡ 1 mod 257 | Full |
| `E_multiple_18_implies_q_mod_19` | If 18 divides E and q ≠ 19, then q ≡ 1 mod 19 | Full |
| `E_multiple_4_implies_c_multiple_4` | If 4 divides E, then 4 divides c | Full |
| `mod_5_implies_c_plus_rho` | If q ≢ 1 mod 5, then c + ρ ≡ 2 mod 4 | Full |
| `height_Q` | Q ≤ 2^(E/2) + 1 | Full |
| `height_c` | c ≤ E/2 | Full |

## Compute reciprocity results

```python
# Include two known even-ρ examples for comparison. Both fail earlier conditions.
selected = sorted(set(keep) | {1093, 3511})
summary = reciprocity_summary(
    selected, all_levels=True, workers=4, silent=True,
    csv_path="q-reciprocity.csv",
)
```

The summary records r, E, ρ, c, the residue for `h = M/E`, and the test values. For odd ρ, `h_candidate` is blank and the report labels the residue as formal. By default the summary computes the grouped result at the largest power of each prime in r. `all_levels=True` includes smaller powers. Levels that require an uncomputed local term or exceed `max_m` (default 128) are marked in the `status` column. Use `blocks=True` to compare the separate cyclotomic factors. That takes more work than the grouped modular calculation.

## Build and sort the report

```python
q_comparison_report(
    "q-conditions.csv", "q-reciprocity.csv", "q-comparison.md",
    qs=selected, normalized_csv="q-normalized.csv",
    sort_by="log_bound_ratio",
)
```

The Markdown has a section per q, a linked table of all selected primes, and a normalized comparison. The normalized Markdown table and CSV sort in ascending order, with blanks at the end. Available keys are `"q"` (default), `"log_bound_ratio"`, `"h_over_Q"`, and `"h_over_q_squared"`. Omit `qs` to include every q in either input file.

For even ρ, start with the logarithmic bound ratio. The candidate needs a ratio below 1. The gap lift checks a further congruence beyond the one used to reconstruct h. When comparing symbols at the same level m, their orders are easier to compare than the raw values of `T mod q`. The per-q table separates the exponent coming from a from the cyclotomic correction.

## Inspect one example

```python
report = explore_reciprocity(433, 2, factor_ideals=True)
gap_candidate(1093)
```

The first command compares the grouped and separate products for s = 2, then checks the actual orders at their prime factors. The second computes the only possible multiplier under the remaining-case hypotheses and reports its failed checks.

## Stream matches or sample condition patterns

```python
pattern = "72817,1,1,1,1,1,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1"
matches = iter_q_matches(
    islice(Primes(), 10^8), example_row=pattern, csv_path="q-matches.csv"
)
for q, flags in matches:
    pass  # Consume the iterator to write the matches.
```

The first number in the example row identifies its source and is ignored
when matching. Use `required={...}` instead of `example_row` for a partial
mask. Cheap checks reject primes before factoring q−1. The scan is
serial and holds no list of candidates. A large input still requires
generating and checking the primes.

```python
selected, types = sample_q_types(
    "q-conditions.csv", limit=10^6, top=5, per_type=2,
    required=dict.fromkeys(Q_CHEAP_CONDITIONS[:6], True),
)
for entry in types:
    print(entry)
```

This consumes the first million data rows, groups identical flag
patterns, and returns the union of the first two examples per selected
pattern. Patterns rank by fewest failed flags, with frequency breaking
ties. Use `ignore=("rho_even",)` to remove a redundant flag from the
ranking while preserving it in the reported pattern. Pass `selected`
to `reciprocity_summary` for the next comparison.

## Native scan through the first 100 million primes

From the repository root:

```sh
make scan_q_native
./build/bin/scan_q_native 2038074743 12 > q-native-100m.csv
```

The numeric limit is inclusive. The scanner prints the prime count and
filter totals to stderr. It first requires the six structural flags.
It writes a row if both order conditions hold or if ρ is at least 2.
The columns match the Sage conditions CSV. Rows arrive in worker order.
Each worker streams its matches and keeps no candidate vector.

For a larger repeated-factor search, append `candidate` to the
command. This applies the first six flags and tests q² divisibility
before factoring q−1. It writes only repeated-factor hits. Append
`wieferich` instead to check all odd primes for repeated factors
without the first six flags. Both modes compute the exact ρ on a hit.
For q above the 32-bit range the square-modulus calculation uses GMP.

```sh
./build/bin/scan_q_native 100000000000 12 candidate > q-native-candidates-1e11.csv
./build/bin/scan_q_native 10000 1 wieferich > q-native-rho-controls.csv
```

The optional fourth argument sets the inclusive lower bound. To extend
the completed scan through q = 10¹¹ without rechecking it:

```sh
./build/bin/scan_q_native 1000000000000 12 candidate 100000000001 > q-native-candidates-1e11-to-1e12.csv
```

## Scan the large-base exponent pairs

Use the native build command above. Run from the repository root.

```sh
./build/bin/scan_q_native prefix 7 3300 3400 --compare > candidates.csv 2> comparison.log
./build/bin/scan_q_native prefix 7 3321924 3321929 > million-digit-band.csv
./build/bin/scan_q_native factor-candidates < candidates.csv > mersenne-factors.csv
```

Each run fixes the odd exponent n₂ ≥ 7. The two endpoints are inclusive
values of j = floor(m₁/n₂). A forced odd q at position j has j+1 bits.
The second command samples a narrow band near one million decimal digits.
Leave stderr on the terminal to see the progress bar, or use
`tail -f comparison.log` in another terminal.

For m₁ = j n₂ + s, the scan takes one exact GMP integer root per reduced
residue s to read the bits of 2^(s/n₂). It extends the prefix if the
last run crosses the computed boundary. The default mode first looks for
a single zero that excludes the whole band. A reported exclusion_zero
is a fractional bit position, numbered from one. Otherwise the scan
checks the individual runs. The required length is (N_min−1)j, where
N_min is the least n₂−e with e positive and even, 2e < n₂, n₂−e ≥ 5,
and gcd(n₂−e,6) = 1.

With `--compare` the scan counts all odd forced bases with even m₁ and
reports the longest run and its position within each mod 240 branch per
s. Ties use the largest j. A count of zero means the branch had no
positions. The passing branch requires M divisible by four and q ≡ 1
modulo 240, or M ≡ 2 modulo four and q ≡ 97 or 193 modulo 240.

Both modes write a candidate CSV row for every qualifying run, including
controls that fail mod 240. The last three flags record the residue test,
q² dividing 2^M−1, and an exact equality of the two power sums with
admissible exponents. The exact check is run only on q² survivors.
Primality of p and q requires separate verification.

The factor command uses the candidate's M = m₁−m₂. It factors each
distinct 2^M−1 once, verifies primality of the factors and checks the
product with multiplicities. For even M it factors the two coprime
halves separately. Each completed M is marked on stderr. Complete
factorization of a million-digit Mersenne number may be computationally
out of reach. You can select one exponent with `./build/bin/scan_q_native factor 60`.

Prefix scans and both factor commands default to the detected logical CPU
count (one if detection is unavailable). Use `--threads K` to limit
concurrent jobs and their memory use. Each prefix worker holds its own
root calculation. With n₂ = 7 there are six independent prefixes, so
six workers cover that scan. An even M provides two independent factor
jobs, and additional exponents can run alongside them. CSV rows arrive
in worker completion order. Progress counts completed prefixes.

```sh
./build/bin/scan_q_native prefix 7 3321924 3321929 --compare --threads 6 > candidates.csv
./build/bin/scan_q_native factor-candidates --threads 12 < candidates.csv > mersenne-factors.csv
```

For counts by residue class, run `./build/bin/scan_q_native 10000000000 24 branch`.
The CSV counts primes passing the original cheap filters, then those also
passing the branch condition. For q ≡ 1 mod 240 this requires 4 | c,
where c = v₂(q−1). For q ≡ 97 or 193 it requires 2^(2u/3) ≡ 1 mod q,
where u = (q−1)/2^c. Candidate mode now applies the branch condition
before testing modulo q². Its stderr cheap count includes that condition.

For class 1 modulo 240, `final_pass` also requires
2^((q−1)/30) ≡ 1 modulo q. This combines the order conditions at 3
and 5. The other two classes retain their `branch_pass` counts.
Candidate mode applies this extra test before the square-modulus test.

## Recover the residue of N modulo q−1

From the project root run `pixi run sage`. To reproduce the saved comparison:

```python
load("scripts/inspect_reciprocity.sage")
qs = set()
for name in ("q-native-100m", "q-native-rho-controls", "q-native-candidates-1e11", "q-native-candidates-1e11-to-1e13"):
    with open("data/" + name + ".csv") as stream:
        qs.update(ZZ(row["q"]) for row in csv.DictReader(stream))
rows = n_remainder_summary(sorted(qs), formal=True, orders=True, silent=True,
                           csv_path="data/q-n-remainders.csv")
```

For one calculation use `n_remainder_window(54001, 2167200, 24079)`.
The arguments are q, M, and a strict integer lower bound on N. The
returned upper bound is inclusive. The first possible N above the lower
bound is recorded only if it lies in the window. Successive values
differ by q−1. Other conditions on N can eliminate further values.

Without `formal=True` the summary skips the window for odd ρ. Formal
rows use the least positive h residue for comparison. With `orders=True`
it also computes the multiplicative order of 2 modulo the odd part u
of q−1. The mandatory divisor of h is ord_u(2)/gcd(ord_u(2), E).
The window test itself uses no order computation modulo u.

## Reverse lift and interactive lab

Start `pixi run jupyter lab labs/prime-power.ipynb` from the project
root. Run Setup and then the section for the calculation of interest.
The notebook includes the reverse and forward lifts, prime scans,
comparison reports, reciprocity checks, Mersenne tables and Lua proofs.

In Sage, after loading `scripts/inspect_reciprocity.sage`:

```python
rows = reverse_n_experiment(reverse_n_branch_primes(10^7), formal=True,
                            csv_path="reverse-n.csv")
trace = reverse_n_lift(54001, 2167200, 1, 24079, 137856, trace=True)
```

The trace stops at the first unique representative or empty interval.
The real-window test checks its representative separately at the
reported precision. first_failed_bit is numbered from 1, so a value
of 9 denotes failure modulo 512. Formal odd-ρ rows are controls.
