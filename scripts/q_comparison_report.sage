# In Sage:
#   load("scripts/q_comparison_report.sage")
#   q_comparison_report("q-conditions.csv", "q-reciprocity.csv", "comparison.md")

from collections import defaultdict
from pathlib import Path
import csv


def _comparison_read_csv(path, required):
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain these columns: {', '.join(required)}")
        return list(reader)


def _comparison_index(rows, key):
    indexed = {}
    for row in rows:
        k = key(row)
        if k in indexed and indexed[k] != row:
            raise ValueError(f"Conflicting CSV rows for {k}.")
        indexed[k] = row
    return indexed


def _comparison_table(out, headers, rows):
    for row in (headers, ["---"] * len(headers), *rows):
        cells = [str(cell).replace("|", r"\|") if cell not in (None, "") else "—"
                 for cell in row]
        out.write("| " + " | ".join(cells) + " |\n")


def _comparison_flag(value):
    if value.lower() in ("1", "true"):
        return True
    if value.lower() in ("0", "false"):
        return False
    raise ValueError(f"Expected a Boolean condition in the CSV, got {value!r}.")


def _comparison_symbols(q, data, row):
    """Express the a-part and correction in the same root B at q."""
    m = ZZ(row["m"])
    F = Zmod(q)
    B = F(2)^(data["E"]//m)
    a_symbol = F(data["a_mod_q"])^((q-1)//m)
    total = F(ZZ(row["T_mod_q"]))
    ka = discrete_log(a_symbol, B, ord=m)
    kt = discrete_log(total, B, ord=m)
    return ka, (kt-ka) % m, kt, m // gcd(m, kt)


def _comparison_ratios(q, data):
    """Return the index and the three ratios shown in the normalized table."""
    if "r" not in data:
        return [None]*4
    R = RealField(80)
    position = R(data["h_residue"])/data["h_modulus"]
    index = (q-1)//data["r"]
    if "h_candidate" not in data:
        return [index, position, None, None]
    h = data["h_candidate"]
    M = h*data["E"]
    log_q = R(q).log()
    height = max((3*M/log_q).log() + R("0.38"), R(10))
    # Left / right in the logarithmic bound on M. Values >=1 fail it.
    lhs = M*R(2).log()/2
    rhs = R("25.2")*log_q^2*height^2 + (R(27)/26).log()
    return [index, position, R(h)/q^2, lhs/rhs]


def _comparison_gap_lift(q, data):
    """Check the first gap congruence beyond the one used to recover h."""
    if "h_candidate" not in data:
        return None, None, None
    h, rho, Q = (data[key] for key in ("h_candidate", "rho", "h_modulus"))
    M = h*data["E"]
    digits = rho+1 if h == 1 else 2*rho
    modulus = q^digits
    A_mod = (power_mod(2, M, Q*modulus)-1)//Q
    return M, digits, (1+power_mod(2, data["c"], modulus)*A_mod) % modulus


def _comparison_prime(out, q, flags, tests):
    """Write one q section and return the two corresponding summary rows."""
    out.write(f"### q {q}\n\n")
    out.write("[All q](#all-q) · [Normalized comparison](#normalized-comparison)\n\n")
    data = {}
    for key in ("r", "E", "rho", "c", "a_mod_q", "h_residue", "h_modulus", "h_candidate"):
        values = {ZZ(row[key]) for row in tests if row.get(key)}
        if len(values) > 1:
            raise ValueError(f"Conflicting values for {key} at q={q}.")
        if values:
            data[key] = values.pop()
    status = "candidate" if "h_candidate" in data else (
        "formal only" if "h_residue" in data else "absent")
    # Since h*2^c*a ≡ -1 mod Q, the saved h also determines a mod Q.
    a_mod_Q = None
    if "h_residue" in data:
        Q = data["h_modulus"]
        a_mod_Q = (-inverse_mod(power_mod(2, data["c"], Q)*data["h_residue"], Q)) % Q
    _comparison_table(out, ("r", "E", "ρ", "c", "a mod Q", "h residue", "Q", "status of h"),
                      [[data.get(k) for k in ("r", "E", "rho", "c")] +
                       [a_mod_Q, data.get("h_residue"), data.get("h_modulus"), status]])
    out.write("\n")
    M, lift_digits, lift_residue = _comparison_gap_lift(q, data)
    lift_status = None
    if M is not None:
        lift_status = "pass" if lift_residue == 0 else "fail"
        _comparison_table(out, ("Candidate M = hE", "Lift modulus", "Residue of 1 + 2^c A", "Gap lift"),
                          [[M, f"q^{lift_digits}", lift_residue, lift_status]])
        out.write("\n")
    if flags:
        _comparison_table(out, ("Failed q conditions", "Passed q conditions"),
                          [[", ".join(k for k in flags if not flags[k]) or "none",
                            ", ".join(k for k in flags if flags[k]) or "none"]])
        out.write("\n")

    test_rows, depth = [], []
    computed = nontrivial = block_levels = block_hits = 0
    for row in tests:
        ka = kc = kt = order = None
        if row["status"] == "computed":
            ka, kc, kt, order = _comparison_symbols(q, data, row)
            computed += 1
            nontrivial += order > 1
            s, m = ZZ(row["s"]), ZZ(row["m"])
            depth.append(f"{m}: {order.valuation(s)}/{m.valuation(s)}")
        if row.get("block_cancellation", "") != "":
            block_levels += 1
            block_hits += bool(row.get("nontrivial_blocks"))
        test_rows.append([row["s"], row["m"], row["status"], row["T_mod_q"],
                          ka, kc, kt, order, row.get("nontrivial_blocks"),
                          row.get("block_cancellation")])
    if tests:
        _comparison_table(out, ("s", "m", "Test", "T mod q", "k(a)", "k(correction)",
                                "k(T)", "order of T", "Nontrivial block indices t",
                                "Block cancellation"), test_rows)
        out.write("\n")

    listed = sum(bool(row["m"]) for row in tests)
    coverage = f"{computed}/{listed}" if listed else None
    link = f"[{q}](#q-{q})"
    all_row = [link, data.get("r"), data.get("E"), data.get("rho"), data.get("c"),
               data.get("a_mod_q"), data.get("h_residue"), data.get("h_modulus"), status,
               M, lift_status,
               sum(not v for v in flags.values()) if flags else None,
               f"{nontrivial}/{computed}" if computed else None, coverage,
               f"{block_hits}/{block_levels}" if block_levels else None]
    normalized_row = [q, status, *_comparison_ratios(q, data), lift_status,
                      coverage, "; ".join(depth) or None]
    return all_row, normalized_row


def q_comparison_report(conditions_csv, reciprocity_csv, output="q-comparison.md",
                        qs=None, normalized_csv=None, sort_by="q"):
    """Join the Sage CSV exports by q and write a Markdown comparison.

    By default includes every q in either file. Use qs=... to select a
    finite subset, or use conditions_csv=None for a reciprocity-only report.
    sort_by orders the normalized table and its optional CSV in ascending
    order, with blanks last. Options: q, h_over_Q, h_over_q_squared,
    log_bound_ratio.
    A skipped reciprocity row stays in the report with its recorded reason.
    """
    sort_columns = {"q": 0, "h_over_Q": 3, "h_over_q_squared": 4, "log_bound_ratio": 5}
    if sort_by not in sort_columns:
        raise ValueError(f"sort_by must be one of: {', '.join(sort_columns)}")
    sort_column = sort_columns[sort_by]
    inputs = {Path(path).resolve() for path in (conditions_csv, reciprocity_csv) if path is not None}
    targets = [Path(path).resolve() for path in (output, normalized_csv) if path is not None]
    if len(set(targets)) != len(targets) or inputs.intersection(targets):
        raise ValueError("Choose distinct paths for the inputs and report outputs.")

    conditions = {}
    if conditions_csv is not None:
        rows = _comparison_read_csv(conditions_csv, ("q",))
        conditions = _comparison_index(rows, lambda row: ZZ(row["q"]))
    rows = _comparison_read_csv(reciprocity_csv, ("q", "s", "m", "status", "T_mod_q"))
    levels = _comparison_index(rows, lambda row: (ZZ(row["q"]), row["s"], row["m"]))
    by_q = defaultdict(list)
    for (q, _, _), row in levels.items():
        by_q[q].append(row)
    primes = set(conditions) | set(by_q)
    if qs is not None:
        primes &= {ZZ(q) for q in qs}

    all_rows, normalized_rows = [], []
    normalized_header = (
        "q", "status of h", "(q−1)/r", "h residue / Q", "h / q²", "log bound ratio",
        "gap lift", "computed / listed", "symbol depth by m",
    )
    with open(output, "w", encoding="utf-8") as out:
        out.write("# Comparing q, h, and the reciprocity tests\n\n")
        out.write("[Individual primes](#individual-primes) · [All q](#all-q) · "
                  "[Normalized comparison](#normalized-comparison)\n\n")
        out.write("## Individual primes\n\n")
        for q in sorted(primes):
            flags = {name: _comparison_flag(value) for name, value in conditions.get(q, {}).items()
                     if name != "q"}
            tests = sorted(by_q.get(q, []), key=lambda row: (ZZ(row["s"] or 0), ZZ(row["m"] or 0)))
            all_row, normalized_row = _comparison_prime(out, q, flags, tests)
            all_rows.append(all_row)
            normalized_rows.append(normalized_row)

        out.write("## All q\n\n")
        _comparison_table(out, ("q", "r", "E", "ρ", "c", "a mod q", "h residue", "Q", "status of h",
                                "Candidate M", "Gap lift", "Failed q flags",
                                "Nontrivial / computed", "Computed / listed",
                                "Levels with block hits / compared"), all_rows)
        out.write("\n## Normalized comparison\n\n")
        normalized_rows.sort(key=lambda row: (row[sort_column] is None,
                                              row[sort_column] or 0, row[0]))
        # Sort the numeric ratios before rounding them for display.
        for row in normalized_rows:
            for index in (3, 4, 5):
                if row[index] is not None:
                    row[index] = row[index].str(digits=5)
        _comparison_table(out, normalized_header,
                          [[f"[{row[0]}](#q-{row[0]})", *row[1:]] for row in normalized_rows])
        out.write(
            "\n**h residue / Q** locates the least positive residue in its full range. "
            "**h / q²** and the **log bound ratio** apply to the candidate branch. "
            "Both ratios must be strictly below 1. The latter is the left side divided "
            "by the right side of the logarithmic bound for M = hE. Values are rounded. "
            "**Gap lift** checks the necessary congruence above the modulus used to recover h.\n\n"
            "For each m = s^b, **k(a)** is the exponent of a^((q−1)/m) in the root "
            "B = 2^(E/m) modulo q. **k(correction)** is the exponent of the remaining "
            "factor, so their sum modulo m gives **k(T)**. In the normalized table, "
            "**m: v/b** gives v = v_s(order(T)) relative to b. Only computed tests enter "
            "that column. Coverage counts the levels listed in the input CSV.\n\n"
            "A nontrivial grouped or separate product certifies an order obstruction "
            "somewhere among its prime factors. A product of 1 can result from "
            "cancellation or from a symbol of insufficient depth. The report uses "
            "the saved residues and needs no factorization of a.\n"
        )
    if normalized_csv is not None:
        with open(normalized_csv, "w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(normalized_header)
            writer.writerows(normalized_rows)
    return Path(output)
