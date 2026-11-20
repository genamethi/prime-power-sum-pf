# Load in Sage with load("scripts/inspect_reciprocity.sage").
# explore_reciprocity(q, s, factor_ideals=True) compares all blocks for s.
# Without factor_ideals it computes the products using only the primes at q.
# gap_candidate(q) reconstructs the only possible M in the even-rho branch.

import csv
import sys
from contextlib import nullcontext
from itertools import islice
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from multiprocessing import get_context


# The column order is also the order of the Boolean entries returned to Sage.
Q_CHEAP_CONDITIONS = (
    "q_gt_3", "mod_3", "mod_8", "mod_240", "c_even_ge_4",
    "r_divides_q_minus_1_over_6", "rho_ge_2", "rho_even",
)
_FORCED_E_CONGRUENCES = ((3, 7), (4, 5), (8, 17), (12, 13), (16, 257), (18, 19))
Q_CONDITIONS = Q_CHEAP_CONDITIONS + (
    "E_rad_r_divides_q_minus_1", "ordinary_divisors_q_minus_1",
    "E_multiple_3_implies_q_mod_9",
) + tuple(f"E_multiple_{d}_implies_q_mod_{ell}" for d, ell in _FORCED_E_CONGRUENCES) + (
    "E_multiple_4_implies_c_multiple_4", "mod_5_implies_c_plus_rho",
    "height_Q", "height_c",
)


def _prime_list(qs):
    """Accept a finite iterable, or one Sage or Python integer."""
    return [ZZ(qs)] if isinstance(qs, (int, Integer)) else [ZZ(q) for q in qs]


def _check_prime(q):
    if not q.is_prime():
        raise ValueError(f"Expected a prime, got {q}.")


def _map_primes(job, qs, workers, batch_size):
    workers, batch_size = int(workers), int(batch_size)
    if min(workers, batch_size) < 1:
        raise ValueError("workers and batch_size must be positive.")
    if workers == 1 or len(qs) < 2:
        return [job(q) for q in qs]
    with ProcessPoolExecutor(max_workers=min(workers, len(qs)),
                             mp_context=get_context("fork")) as pool:
        return list(pool.map(job, qs, chunksize=batch_size))


def _write_rows(header, rows, silent, csv_path):
    if not silent:
        writer = csv.writer(sys.stdout, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    if csv_path is not None:
        with open(csv_path, "w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)


def _rho_from_multiple(q, multiple):
    """Compute v_q(2^multiple - 1) by modular exponentiation."""
    rho, modulus = 0, q
    while power_mod(2, multiple, modulus) == 1:
        rho += 1
        modulus *= q
    return rho


def _order_data(q, factors):
    """Strip prime factors from the known order multiple q-1."""
    r = q - 1
    for ell, _ in factors:
        while r % ell == 0 and power_mod(2, r // ell, q) == 1:
            r //= ell
    return r, lcm(2, r)


def _q_basic_flags(q):
    c = (q-1).valuation(2)
    return [q > 3, q % 3 == 1, q % 8 == 1, q % 240 in (1, 97, 193),
            c >= 4 and c % 2 == 0]


def _q_condition_row(q, cheap, cheap_row=None):
    #_check_prime(q)
    if q == 2:
        # The order of 2 modulo 2 is undefined.
        return [False] * len(Q_CHEAP_CONDITIONS if cheap else Q_CONDITIONS)
    c = (q - 1).valuation(2)
    # (q-1)/r is prime to q, so using q-1 gives the same valuation as r.
    rho = _rho_from_multiple(q, q - 1) if cheap_row is None else None
    row = cheap_row if cheap_row is not None else _q_basic_flags(q) + [
        (q-1) % 6 == 0 and power_mod(2, (q-1)//6, q) == 1,
        rho >= 2, rho % 2 == 0,
    ]
    if cheap:
        return row

    if rho is None:
        rho = 1 if not row[6] else _rho_from_multiple(q, q-1)
    factors = factor(q - 1)
    r, E = _order_data(q, factors)
    rad_r = prod(ell for ell, _ in factors if r % ell == 0)
    ordinary = [ell for ell, _ in factors if ell != 2
                and power_mod(2, ell-1, ell^2) != 1]
    Q = q^rho
    row += [(q-1) % (E*rad_r) == 0,
            all((q-1)//ell % r == 0 for ell in ordinary),
            E % 3 != 0 or q % 9 == 1]
    row += [E % d != 0 or q == ell or q % ell == 1 for d, ell in _FORCED_E_CONGRUENCES]
    # ceil(log_2(Q-1)) = bit length of Q-2.
    row += [E % 4 != 0 or c % 4 == 0,
            q % 5 == 1 or (c+rho) % 4 == 2,
            E >= 2*(Q-2).nbits(), E >= 2*c]
    return row


def q_conditions(qs, silent=False, cheap=False, csv_path=None, workers=1, batch_size=64):
    """Return Boolean rows, one per input prime in its original order.

    Use Q_CONDITIONS (or Q_CHEAP_CONDITIONS with cheap=True) for column names.
    Implications pass when their hypotheses are false. The forced prime
    congruences also require q distinct from the prime in the modulus.
    height_Q tests q^rho <= 2^(E/2)+1. height_c tests c <= E/2.
    Use a finite iterable, for example islice(Primes(), 1000).
    Every input produces a row. silent=True suppresses the printed CSV
    table, but csv_path still writes a file if given.
    """
    qs = _prime_list(qs)
    names = Q_CHEAP_CONDITIONS if cheap else Q_CONDITIONS
    rows = _map_primes(partial(_q_condition_row, cheap=cheap), qs, workers, batch_size)
    if not silent or csv_path is not None:
        _write_rows(("q",) + names,
                    [[q] + [int(flag) for flag in row] for q, row in zip(qs, rows)],
                    silent, csv_path)
    return rows


def iter_q_matches(qs, required=None, example_row=None, csv_path=None):
    """Yield (q, row) from a prime iterator and optionally write only matches."""
    required = dict(required or {})
    if example_row is not None:
        bits = [int(cell) for cell in example_row.split(",")]
        if required or len(bits) != len(Q_CONDITIONS)+1 or any(x not in (0, 1) for x in bits[1:]):
            raise ValueError("Supply either required or an example row of q and 21 flags.")
        required = dict(zip(Q_CONDITIONS, bits[1:]))
    indices = [(Q_CONDITIONS.index(name), value) for name, value in required.items()]
    basic = [(i, value) for i, value in indices if i < 5]
    cheap = [(i, value) for i, value in indices if i < len(Q_CHEAP_CONDITIONS)]
    with (open(csv_path, "w", newline="") if csv_path is not None else nullcontext()) as stream:
        writer = csv.writer(stream) if stream is not None else None
        if writer is not None:
            writer.writerow(("q",) + Q_CONDITIONS)
        for q in qs:
            q = ZZ(q)
            if basic and q > 2:
                flags = _q_basic_flags(q)
                if any(flags[i] != value for i, value in basic):
                    continue
            row = _q_condition_row(q, True)
            if any(row[i] != value for i, value in cheap):
                continue
            row = _q_condition_row(q, False, cheap_row=row)
            if any(row[i] != value for i, value in indices):
                continue
            if writer is not None:
                writer.writerow([q] + [int(flag) for flag in row])
            yield q, row


def sample_q_types(csv_path, top=5, per_type=5, limit=None, required=None, ignore=()):
    """Return (union of sampled primes, pattern summaries), ranked by failures."""
    groups = {}
    required = dict(required or {})
    with open(csv_path, newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        if not header or header[0] != "q":
            raise ValueError("Expected a q-conditions CSV.")
        names = header[1:]
        mask = [(names.index(name), str(int(value))) for name, value in required.items()]
        ranked = [i for i, name in enumerate(names) if name not in ignore]
        for q, *cells in islice(reader, limit):
            if any(cells[i] != value for i, value in mask):
                continue
            pattern = tuple(cells)
            count, qs = groups.setdefault(pattern, [0, []])
            groups[pattern][0] = count+1
            if len(qs) < per_type:
                qs.append(ZZ(q))
    patterns = sorted(groups, key=lambda p: (
        sum(p[i] == "0" for i in ranked), -groups[p][0], p))[:top]
    types = [dict(count=groups[p][0], qs=groups[p][1],
                  failed=[name for name, cell in zip(names, p) if cell == "0"])
             for p in patterns]
    return sorted({q for entry in types for q in entry["qs"]}), types


def cyclotomic_block(z, E, m, t):
    """Relative norm with the choice of root satisfying zeta_t^(E/m) = z."""
    h = t // m
    root = z^inverse_mod(E // t, m)
    # Inclusion-exclusion over divisors coprime to m.
    return prod(
        (2^(h // d) - root^inverse_mod(d, m))^moebius(d)
        for d in divisors(h) if gcd(d, m) == 1 and moebius(d) != 0
    )


def inspect_reciprocity(q, m, t=None, factor_ideals=False):
    """Compare the local contributions at q with the selected ideal product."""
    q, m = ZZ(q), ZZ(m)
    if q <= 2 or not q.is_prime() or m < 2 or len(factor(m)) != 1:
        raise ValueError("Use an odd prime q and a prime power m.")
    s = factor(m)[0][0]
    r = Zmod(q)(2).multiplicative_order()
    E = lcm(2, r)
    if r % m:
        raise ValueError("m must divide ord_q(2).")
    if q % (2*m if s == 2 else s*m) != 1:
        raise ValueError("The local term at s needs a separate calculation.")

    K = CyclotomicField(m)
    z = K.gen()
    B = 2^(E // m)
    element = B - z
    if t is not None:
        t = ZZ(t)
        if t <= 0 or E % t or (E // t) % s == 0:
            raise ValueError("t must divide E and contain its full s-part.")
        element = cyclotomic_block(z, E, m, t)

    at_q = []
    for i in range(1, m):
        if gcd(i, m) != 1:
            continue
        P = K.ideal(q, z - power_mod(B, i, q))
        F = P.residue_field()
        v = element.valuation(P)
        unit = F(element / q^v)
        residue = ((-1)^v * unit)^((q - 1) // m)
        k = discrete_log(residue, F(z), ord=m)
        at_q.append(dict(i=i, valuation=v, unit=ZZ(unit), exponent=k))

    exponent = sum(row["exponent"] for row in at_q) % m
    result = dict(q=q, r=r, E=E, m=m, t=t, at_q=at_q,
                  exponent=exponent, T=power_mod(B, exponent, q))
    if factor_ideals:
        factors = []
        for P, v in K.ideal(element).factor():
            ell = P.smallest_integer()
            if ell in (q, s):
                continue
            F = P.residue_field()
            residue = F(q)^((F.cardinality() - 1) // m)
            k = discrete_log(residue, F(z), ord=m)
            order = Zmod(ell)(q).multiplicative_order()
            order_power = Zmod(ell^v)(q).multiplicative_order()
            factors.append(dict(ell=ell, valuation=v, exponent=k,
                                order=order, order_power=order_power))
        if sum(row["valuation"] * row["exponent"] for row in factors) % m != exponent:
            raise ArithmeticError("The product at the prime factors differs from the product at q.")
        result["factors"] = factors
    return result


def explore_reciprocity(q, s, b=None, factor_ideals=False):
    """Print a comparison of all cyclotomic factors for one prime s | ord_q(2).

    By default m is the full power of s in ord_q(2). Set b for a smaller
    level. With factor_ideals=True, show which primes witness s | ord_ell(q)
    and which witnesses are missed by the m-th power-residue symbol.
    """
    q, s = ZZ(q), ZZ(s)
    if not s.is_prime():
        raise ValueError("s must be prime.")
    if b is None:
        b = Zmod(q)(2).multiplicative_order().valuation(s)
    if b < 1:
        raise ValueError("s must divide ord_q(2).")
    m = s^b
    grouped = inspect_reciprocity(q, m)
    E = grouped["E"]
    indices = [t for t in divisors(E) if (E // t) % s != 0]
    blocks = [inspect_reciprocity(q, m, t, factor_ideals) for t in indices]
    if sum(row["exponent"] for row in blocks) % m != grouped["exponent"]:
        raise ArithmeticError("The block products do not multiply to the grouped product.")

    print(f"q = {q}, r = {grouped['r']}, E = {E}, m = {m}")
    print("k means product = zeta_m^k, with zeta_m reducing to 2^(E/m) at q.")
    print("factor            k     T mod q     conclusion")
    for label, row in [("grouped", grouped)] + [(f"Phi_{b['t']}(2)", b) for b in blocks]:
        k = row["exponent"]
        conclusion = f"some ell has {s} | ord_ell(q)" if k else "product is 1"
        print(f"{label:<17} {str(k):>3} {str(row['T']):>11}     {conclusion}")
    if grouped["exponent"] == 0 and any(row["exponent"] for row in blocks):
        print("The nontrivial block products cancel in the grouped product.")

    if factor_ideals:
        print("\nIndividual factors. The depth column is v_s(ell-1).")
        print("t       ell^f                       ord_ell(q)    depth     k   conclusion")
        for block in blocks:
            for row in block["factors"]:
                ell, order, k = row["ell"], row["order"], row["exponent"]
                depth = (ell-1).valuation(s)
                if k:
                    conclusion = "detected by symbol"
                elif order % s == 0:
                    min_b = depth - order.valuation(s) + 1
                    conclusion = f"order divisible by s, needs b >= {min_b} at this ell"
                elif gcd(row["order_power"], E) > 1:
                    conclusion = f"other order obstruction: gcd(ord_ell^f(q), E) = {gcd(row['order_power'], E)}"
                else:
                    conclusion = "no order obstruction"
                label = f"{ell}^{row['valuation']}"
                print(f"{str(block['t']):<7} {label:<24} {str(order):>13} "
                      f"{str(depth):>8} {str(k):>5}   {conclusion}")
    return dict(grouped=grouped, blocks=blocks)


def _h_residue(q, E, rho, c):
    Q = q^rho
    a = (power_mod(2, E, Q^2) - 1) // Q
    return (-inverse_mod(power_mod(2, c, Q) * a, Q)) % Q, Q, a


RECIPROCITY_COLUMNS = (
    "q", "r", "E", "rho", "c", "a_mod_q", "h_residue", "h_modulus",
    "h_candidate", "s", "m", "status", "T_mod_q", "nontrivial_blocks", "block_cancellation",
)


def _reciprocity_rows(q, max_m, all_levels, blocks):
    #_check_prime(q)
    if q == 2:
        return [dict(q=q, status="order undefined")]
    factors = factor(q-1)
    r, E = _order_data(q, factors)
    rho = _rho_from_multiple(q, q-1)
    c = (q-1).valuation(2)
    h, Q, a = _h_residue(q, E, rho, c)
    base = dict(q=q, r=r, E=E, rho=rho, c=c, a_mod_q=a % q,
                h_residue=h, h_modulus=Q, h_candidate=h if q > 3 and rho % 2 == 0 else None)
    rows = []
    for s, _ in factors:
        full_depth = r.valuation(s)
        if full_depth == 0:
            continue
        levels = range(1, full_depth+1) if all_levels else (full_depth,)
        for b in levels:
            m = s^b
            row = dict(base, s=s, m=m)
            rows.append(row)
            if q % (2*m if s == 2 else s*m) != 1:
                row["status"] = "local term at s not covered"
                continue
            if max_m is not None and m > max_m:
                row["status"] = "above max_m"
                continue
            F = Zmod(q)
            B = F(2)^(E//m)
            # (B^(m/s)-1) Phi'_m(B) = m/B. The sign is needed for odd rho.
            W = F((-1)^rho * a)*B/F(m)
            for i in range(2, m):
                if i % s:
                    W *= (B-B^i)^inverse_mod(i, m)
            value = ZZ(W^((q-1)//m))
            row.update(status="computed", T_mod_q=value)
            if blocks:
                nontrivial = []
                product = F(1)
                for t in divisors(E):
                    if (E//t) % s:
                        block = inspect_reciprocity(q, m, t)
                        product *= F(block["T"])
                        if block["exponent"]:
                            nontrivial.append(t)
                if product != F(value):
                    raise ArithmeticError("The block products do not match the grouped product.")
                row.update(nontrivial_blocks=" ".join(map(str, nontrivial)),
                           block_cancellation=value == 1 and bool(nontrivial))
    return rows or [dict(base, status="no prime power in r")]


def reciprocity_summary(qs, silent=False, csv_path=None, workers=1, batch_size=8,
                        max_m=128, all_levels=False, blocks=False):
    """Return dictionaries keyed by RECIPROCITY_COLUMNS (optional fields may be absent).

    Use all_levels=True to include smaller prime-power levels of r. The
    grouped calculation uses modular arithmetic. Set blocks=True to compare
    separate relative norms, without factoring their ideal products.
    h_residue is the formal residue for every odd q, while h_candidate
    is populated only for q>3 and even rho. Such a candidate can still fail
    other restrictions. Pass max_m=None to remove the work cap.
    """
    qs = _prime_list(qs)
    job = partial(_reciprocity_rows, max_m=max_m, all_levels=all_levels, blocks=blocks)
    rows = [row for group in _map_primes(job, qs, workers, batch_size) for row in group]
    _write_rows(RECIPROCITY_COLUMNS,
                [[row.get(key) for key in RECIPROCITY_COLUMNS] for row in rows]
                if not silent or csv_path is not None else (),
                silent, csv_path)
    return rows


def gap_candidate(q, exact_bits=0):
    """Recover h=M/E modulo q^rho and check one higher congruence.

    Results are conditional on the remaining-case hypotheses of the notes.
    With exact_bits > 0, also compare the gap equation exactly when M fits
    the given bit budget. No search over M or N is needed.
    """
    q = ZZ(q)
    if q <= 3 or not q.is_prime():
        raise ValueError("Use a prime q > 3.")
    r = Zmod(q)(2).multiplicative_order()
    E = lcm(2, r)
    c = (q - 1).valuation(2)
    rho = _rho_from_multiple(q, E)
    result = dict(q=q, r=r, E=E, rho=rho, m2=c)
    if rho % 2:
        result["status"] = "odd rho excludes the remaining branch"
        return result
    h, Q, a = _h_residue(q, E, rho, c)
    result.update(h=h, M=h*E, a_mod_Q=a)
    failures = []
    if q % 240 not in (1, 97, 193) or c % 2 or (q-1) % (3*r):
        failures.append("earlier congruences on q and r")
    if h >= q^2:
        failures.append("h is outside 1 <= h < q^2")
    M = h*E
    interval = RealIntervalField(96)
    log_q = interval(q).log()
    height = max(((3*M/log_q).log() + QQ(38)/100).upper(), 10)
    lhs = M * interval(2).log()/2
    rhs = (QQ(252)/10) * log_q^2 * interval(height)^2 + interval(QQ(27)/26).log()
    if lhs.lower() >= rhs.upper():
        failures.append("M exceeds the logarithmic bound")
    if M % 4 == 0 and (q % 5 != 1 or c % 4):
        failures.append("the modulo-5 conditions for 4 | M")
    if M % 4 == 2 and (q % 5 not in (2, 3) or (c + rho) % 4 != 2):
        failures.append("the modulo-5 conditions for M = 2 mod 4")
    if r % 2 == 0 and q % 8 == 1 and power_mod(h, (q-1)//2, q) != 1:
        failures.append("h is not a square modulo q")
    if r % 3 == 0 and (q-1) % (3*r) == 0 and power_mod(h, (q-1)//3, q) != 1:
        failures.append("h is not a cube modulo q")
    if r % 4 == 0 and q % 8 == 1 and c % 4 == 0 and power_mod(2*h, (q-1)//4, q) != 1:
        failures.append("2h is not a fourth power modulo q")

    digits = min(2*rho, (2*h-1)*rho + 1)
    modulus = q^digits
    A = (power_mod(2, h*E, modulus*Q) - 1) // Q
    residue = (1 + power_mod(2, c, modulus)*A) % modulus
    result.update(next_modulus=modulus, gap_residue=residue)
    if residue:
        failures.append(f"q^N would not be 0 modulo q^{digits}")
    result["failures"] = failures
    if exact_bits and h*E + c + 1 <= exact_bits:
        gap = 1 + 2^c * ((2^(h*E) - 1) // Q)
        N = gap.valuation(q)
        result["exact_gap"] = f"q^{N}" if gap == q^N else "not a power of q"
        result["N"] = N if gap == q^N else None
    return result
