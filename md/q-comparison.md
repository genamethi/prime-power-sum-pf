# Comparing q, h, and the reciprocity tests

[Individual primes](#individual-primes) · [All q](#all-q) · [Normalized comparison](#normalized-comparison)

## Individual primes

## All q

| q | r | E | ρ | c | a mod q | h residue | Q | status of h | Candidate M | Gap lift | Failed q flags | Nontrivial / computed | Computed / listed | Levels with block hits / compared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Normalized comparison

| q | status of h | (q−1)/r | h residue / Q | h / q² | log bound ratio | gap lift | computed / listed | symbol depth by m |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

**h residue / Q** locates the least positive residue in its full range. **h / q²** and the **log bound ratio** apply to the candidate branch. Both ratios must be strictly below 1. The latter is the left side divided by the right side of the logarithmic bound for M = hE. Values are rounded. **Gap lift** checks the necessary congruence above the modulus used to recover h.

For each m = s^b, **k(a)** is the exponent of a^((q−1)/m) in the root B = 2^(E/m) modulo q. **k(correction)** is the exponent of the remaining factor, so their sum modulo m gives **k(T)**. In the normalized table, **m: v/b** gives v = v_s(order(T)) relative to b. Only computed tests enter that column. Coverage counts the levels listed in the input CSV.

A nontrivial grouped or separate product certifies an order obstruction somewhere among its prime factors. A product of 1 can result from cancellation or from a symbol of insufficient depth. The report uses the saved residues and needs no factorization of a.
