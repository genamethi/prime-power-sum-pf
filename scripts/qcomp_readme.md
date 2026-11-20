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
g++ -O3 -std=c++17 -pthread scripts/scan_q_native.cpp -lflint -lgmpxx -lgmp -o /tmp/scan_q_native
/tmp/scan_q_native 2038074743 12 > q-native-100m.csv
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
/tmp/scan_q_native 100000000000 12 candidate > q-native-candidates-1e11.csv
/tmp/scan_q_native 10000 1 wieferich > q-native-rho-controls.csv
```

The optional fourth argument sets the inclusive lower bound. To extend
the completed scan through q = 10¹¹ without rechecking it:

```sh
/tmp/scan_q_native 1000000000000 12 candidate 100000000001 > q-native-candidates-1e11-to-1e12.csv
```
