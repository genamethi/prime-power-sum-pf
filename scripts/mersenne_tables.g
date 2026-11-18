# GAP 4; no optional packages. Run from this directory:
#   ~/fluid/sauce/gap/gap -A -q -b mersenne_tables.g
# To change the bounds without editing the script, define PPConfig before Read.
# Setting PPLibraryOnly := true loads the functions without running the report.

PPValuation := function(n, prime)
    local exponent;
    if n <= 0 then Error("valuation requires a positive integer"); fi;
    exponent := 0;
    while n mod prime = 0 do
        n := QuoInt(n, prime);
        exponent := exponent + 1;
    od;
    return exponent;
end;

# Compute v_prime(base^exponent - 1) without constructing base^exponent.
PPPowerValuation := function(base, exponent, prime)
    local valuation, modulus;
    valuation := 0;
    modulus := prime;
    while PowerModInt(base, exponent, modulus) = 1 do
        valuation := valuation + 1;
        modulus := modulus * prime;
    od;
    return valuation;
end;

PPProvenPrime := function(n)
    local proof;
    if not IsPrimeInt(n) then return false; fi;
    if n < 2^64 then return true; fi;
    proof := PrimalityProof(n);
    if proof = fail or proof = false then
        Error("primality not certified for ", n);
    fi;
    return PrimalityProof_Verify(n, proof);
end;

PPFactorization := function(n)
    local factors;
    if n = 1 then return []; fi;
    factors := FactorsInt(n);
    if Product(factors) <> n or not ForAll(Set(factors), PPProvenPrime) then
        Error("factorization not certified");
    fi;
    return Collected(factors);
end;

PPFactorText := function(factors)
    if IsEmpty(factors) then return "1"; fi;
    return JoinStringsWithSeparator(List(factors,function(pair)
        if pair[2]=1 then return String(pair[1]); fi;
        return Concatenation(String(pair[1]),"^",String(pair[2]));
    end)," * ");
end;

# Primes ell | q-1 for which the valuation identity in the notes applies.
PPOrdinaryDivisors := function(q)
    return Filtered(Set(FactorsInt(q-1)), ell -> ell <> 2 and
        PPPowerValuation(2, OrderMod(2, ell), ell) = 1);
end;

PPBaseData := function(q)
    local r;
    r := OrderMod(2, q);
    return rec(q := q, c := PPValuation(q-1, 2), r := r,
        E := Lcm(2, r), e := PPPowerValuation(2, r, q),
        ordinary := PPOrdinaryDivisors(q));
end;

# These restrictions depend on q alone, since E divides M.
PPBaseCompatible := function(b)
    local condition, Q, ceiling;
    if b.e < 2 or b.e mod 2 <> 0 then return false; fi;
    if not ForAll(b.ordinary, ell -> QuoInt(b.q-1, ell) mod b.r = 0) then
        return false;
    fi;
    for condition in [[3,7], [4,5], [8,17], [12,13], [16,257], [18,19]] do
        if b.E mod condition[1] = 0 and b.q mod condition[2] <> 1 then
            return false;
        fi;
    od;
    if b.E mod 4 = 0 and b.c mod 4 <> 0 then return false; fi;
    if b.q mod 5 <> 1 and (b.c+b.e) mod 4 <> 2 then return false; fi;
    Q := b.q^b.e;
    ceiling := LogInt(Q-1, 2);
    if 2^ceiling < Q-1 then ceiling := ceiling+1; fi;
    return 2*ceiling <= b.E and 2*b.c <= b.E;
end;

PPScan := function(limit)
    local counts, candidates, residue, q, b;
    counts := [0,0,0,0,0,0];
    candidates := [];
    # Only these residue classes can occur; ranges are inclusive.
    for residue in [1,97,193] do
        if residue > limit then continue; fi;
        for q in [residue,residue+240..residue+240*QuoInt(limit-residue,240)] do
            if q < 97 or not PPProvenPrime(q) then continue; fi;
            counts[1] := counts[1]+1;
            if PPValuation(q-1,2) mod 2 <> 0 then continue; fi;
            counts[2] := counts[2]+1;
            if PowerModInt(2,QuoInt(q-1,6),q) <> 1 then continue; fi;
            counts[3] := counts[3]+1;
            if PowerModInt(2,q-1,q^2) <> 1 then continue; fi;
            counts[4] := counts[4]+1;
            b := PPBaseData(q);
            if b.e mod 2 <> 0 then continue; fi;
            counts[5] := counts[5]+1;
            if not PPBaseCompatible(b) then continue; fi;
            counts[6] := counts[6]+1;
            Add(candidates,b);
        od;
    od;
    SortBy(candidates, b -> b.q);
    return rec(counts := counts, candidates := candidates);
end;

# Structural rows retain the first three q conditions, with a bound on E.
# Exact valuations and subsequent failures are recorded, not used to erase rows.
PPStructuralScan := function(config)
    local counts,rows,residue,q,r,E;
    counts := [0,0,0,0,0];
    rows := [];
    for residue in [1,97,193] do
        if residue > config.qMax then continue; fi;
        for q in [residue,residue+240..residue+240*QuoInt(config.qMax-residue,240)] do
            if q<97 or not PPProvenPrime(q) then continue; fi;
            counts[1] := counts[1]+1;
            if PPValuation(q-1,2) mod 2<>0 then continue; fi;
            counts[2] := counts[2]+1;
            if PowerModInt(2,QuoInt(q-1,6),q)<>1 then continue; fi;
            counts[3] := counts[3]+1;
            if PowerModInt(2,q-1,q^2)=1 then counts[5] := counts[5]+1; fi;
            r := OrderMod(2,q);
            E := Lcm(2,r);
            if E>config.factorE then continue; fi;
            counts[4] := counts[4]+1;
            Add(rows,PPBaseData(q));
        od;
    od;
    SortBy(rows,b -> b.q);
    return rec(counts := counts,rows := rows);
end;

PPFailureText := function(b)
    local failures,ell,condition;
    failures := [];
    if b.e mod 2<>0 then Add(failures,"rho odd"); fi;
    for ell in b.ordinary do
        if QuoInt(b.q-1,ell) mod b.r<>0 then
            Add(failures,Concatenation("r does not divide (q-1)/",String(ell)));
        fi;
    od;
    for condition in [[3,7],[4,5],[8,17],[12,13],[16,257],[18,19]] do
        if b.E mod condition[1]=0 and b.q mod condition[2]<>1 then
            Add(failures,Concatenation("q != 1 mod ",String(condition[2]),
                " (",String(condition[1])," divides E)"));
        fi;
    od;
    if b.E mod 4=0 and b.c mod 4<>0 then Add(failures,"4 divides E but c != 0 mod 4"); fi;
    if b.q^2 > 2^(b.E/2)+1 then Add(failures,"q^2 > 2^(E/2)+1"); fi;
    if 2*b.c > b.E then Add(failures,"c > E/2"); fi;
    if IsEmpty(failures) then return "none in these tests"; fi;
    return JoinStringsWithSeparator(failures,"; ");
end;

PPFactorsOfA := function(b)
    local a, result, pair, ell, t;
    a := QuoInt(2^b.E-1,b.q^b.e);
    if a*b.q^b.e <> 2^b.E-1 then Error("incorrect q valuation"); fi;
    result := [];
    for pair in PPFactorization(a) do
        ell := pair[1];
        t := OrderMod(2,ell);
        Add(result,rec(ell := ell, f := pair[2], t := t,
            j := QuoInt(ell-1,t), orderQ := OrderMod(b.q,ell),
            orderQPower := OrderMod(b.q,ell^pair[2]),
            removalPrimes := Filtered(List(PPFactorization(QuoInt(ell-1,t)),p -> p[1]),
                s -> Gcd(s,6*b.E)=1),
            phiIndices := Filtered(DivisorsInt(b.E), d ->
                Value(CyclotomicPolynomial(Rationals,d),2) mod ell = 0)));
    od;
    return result;
end;

PPSelectedProduct := function(factors,N)
    return Product(Filtered(factors,f -> Gcd(f.j,N)=1), f -> f.ell^f.f);
end;

# Arithmetic Frobenius in Q(zeta_m). Here m=N is odd, or m=t divides ell-1.
PPFrobenius := function(ell,m)
    if Gcd(ell,m)<>1 then return ["ramified","-"]; fi;
    if m=1 then return [1,1]; fi;
    return [ell mod m,OrderMod(ell,m)];
end;

PPLine := function(stream,line)
    WriteLine(stream,line);
end;

PPRow := function(stream,entries)
    PPLine(stream,Concatenation("| ",JoinStringsWithSeparator(
        List(entries,String)," | ")," |"));
end;

PPDetails := function(stream,b,config)
    local factors,N,f,frob,s,primes,nValues,L,D,selected,ordinaryPass;
    PPLine(stream,Concatenation("\n## q = ",String(b.q),"\n"));
    PPRow(stream,["q","c = v₂(q−1)","rho = v_q(2^r−1)","r","E","ord_(q²)(2)"]);
    PPRow(stream,["---","---","---","---","---","---"]);
    PPRow(stream,[b.q,b.c,b.e,b.r,b.E,OrderMod(2,b.q^2)]);
    PPLine(stream,Concatenation("\nFailed additional conditions: ",PPFailureText(b),"."));
    PPLine(stream,Concatenation("\nq−1 = ",PPFactorText(PPFactorization(b.q-1)),
        "; ordinary odd prime divisors of q−1: ",String(b.ordinary),"."));
    Print("Factoring E = ",b.E," for q = ",b.q,"\n");
    factors := PPFactorsOfA(b);
    PPLine(stream,"\n### Factors of a = (2^E − 1)/q^rho\n");
    PPLine(stream,Concatenation("a = ",String(QuoInt(2^b.E-1,b.q^b.e)),".\n"));
    PPRow(stream,["ℓ","v_ℓ(a)","t = ord_ℓ(2)","j = (ℓ−1)/t",
        "ord_ℓ(q)","ord_(ℓ^f)(q)","indices d with ℓ dividing Φ_d(2)","ℓ > q",
        "ℓ divides q−1","S_ℓ: primes s that can remove ℓ^f"]);
    PPRow(stream,List([1..10],x -> "---"));
    for f in factors do
        PPRow(stream,[f.ell,f.f,f.t,f.j,f.orderQ,f.orderQPower,f.phiIndices,
            f.ell>b.q,(b.q-1) mod f.ell=0,f.removalPrimes]);
    od;
    L := Lcm(List(factors,f -> f.orderQPower));
    PPLine(stream,Concatenation("\nL = ord_a(q) = ",String(L)," = ",
        PPFactorText(PPFactorization(L)),". gcd(L,6E) = ",String(Gcd(L,6*b.E)),"."));
    PPLine(stream,"\n### Selection by N\n");
    PPRow(stream,["N","prime factors of N","prime powers selected for D_E","D_E",
        "D_E divides uN","L divides N","ordinary divisors of q−1 avoid N"]);
    PPRow(stream,List([1..7],x -> "---"));
    nValues := Filtered([5..config.nMax],N -> Gcd(N,6*b.E)=1);
    for N in nValues do
        D := PPSelectedProduct(factors,N);
        selected := Filtered(factors,f -> Gcd(f.j,N)=1);
        ordinaryPass := ForAll(b.ordinary,ell -> N mod ell<>0);
        PPRow(stream,[N,Set(FactorsInt(N)),PPFactorText(List(selected,f -> [f.ell,f.f])),
            D,(QuoInt(b.q-1,2^b.c)*N) mod D=0,N mod L=0,ordinaryPass]);
    od;
    PPLine(stream,"\n### Frobenius in Q(ζ_s)\n");
    PPRow(stream,["s","ℓ","Frob: ℓ mod s","residue degree","ramification index",
        "number of prime ideals above ℓ","passes the s condition for D_E"]);
    PPRow(stream,List([1..7],x -> "---"));
    primes := Set(Concatenation(Concatenation(List(nValues,FactorsInt)),
        Concatenation(List(factors,f -> f.removalPrimes))));
    for s in primes do
        for f in factors do
            frob := PPFrobenius(f.ell,s);
            if f.ell=s then
                PPRow(stream,[s,f.ell,frob[1],1,s-1,1,true]);
            else
                PPRow(stream,[s,f.ell,frob[1],frob[2],1,QuoInt(s-1,frob[2]),
                    (f.ell-1) mod s<>0]);
            fi;
        od;
    od;
end;

PPReport := function(config)
    local scan,stream,labels,i,b;
    scan := PPStructuralScan(config);
    stream := OutputTextFile(config.output,false);
    SetPrintFormattingStatus(stream,false);
    PPLine(stream,"# Mersenne factors, orders and Frobenius\n");
    PPLine(stream,Concatenation("GAP ",GAPInfo.Version,". Prime q range: 97 ≤ q ≤ ",
        String(config.qMax)," (inclusive)."));
    PPLine(stream,"The range is a computational choice, not an upper bound from the notes.");
    PPLine(stream,"Structural tables for the remaining case of `prime_power_notes.tex`: A > 1, q > 3, q does not divide M.");
    PPLine(stream,"Retain q mod 240 in {1,97,193}, c=v₂(q−1) even, and r=ord_q(2) dividing (q−1)/6. Bound E=lcm(2,r) for complete factorization. Subsequent failures are displayed, not filtered out.");
    PPLine(stream,"Set rho=v_q(2^r−1), a=(2^E−1)/q^rho, u=(q−1)/2^c. In a solution, rho=n₁ must be even. An odd rho is recorded as a failure, not renamed n₁.\n");
    PPLine(stream,"## Contents\n\n- [Coverage](#counts-after-successive-conditions)\n- [Overview](#overview)\n- [Table definitions](#table-definitions-and-coverage)\n- [Frobenius and the selected factors](#frobenius-and-the-selected-factors)\n- [Reproduction](#reproduction)");
    for b in scan.rows do
        PPLine(stream,Concatenation("- [q = ",String(b.q),"](#q--",String(b.q),")"));
    od;
    PPLine(stream,"\n## Counts after successive conditions\n");
    PPRow(stream,["Condition retained cumulatively","Number of q"]);
    PPRow(stream,["---","---"]);
    labels := ["q prime; q mod 240 in {1,97,193}","c=v₂(q−1) even (hence c≥4)",
        "2^((q−1)/6) = 1 modulo q; equivalently r divides (q−1)/6",
        Concatenation("E ≤ ",String(config.factorE),"; rows tabulated below")];
    for i in [1..Length(labels)] do PPRow(stream,[labels[i],scan.counts[i]]); od;
    PPLine(stream,Concatenation("\nAmong all ",String(scan.counts[3]),
        " q passing the first three conditions (without the E bound), ",
        String(scan.counts[5])," satisfy 2^(q−1)=1 modulo q²."));
    if scan.counts[5]=0 then
        PPLine(stream,"Thus none of the present rows is a candidate for the full system.");
    fi;
    PPLine(stream,"\n## Overview\n");
    PPRow(stream,["q","q mod 240","c","r","E","rho","ord_(q²)(2)","failed additional conditions"]);
    PPRow(stream,List([1..8],x -> "---"));
    for b in scan.rows do
        PPRow(stream,[b.q,b.q mod 240,b.c,b.r,b.E,b.e,OrderMod(2,b.q^2),PPFailureText(b)]);
    od;
    for b in scan.rows do PPDetails(stream,b,config); od;
    PPLine(stream,"\n## Table definitions and coverage\n");
    PPLine(stream,Concatenation("The q rows exhaust the retained conditions with q ≤ ",String(config.qMax),
        " and E ≤ ",String(config.factorE),". The N tables cover every 5 ≤ N ≤ ",
        String(config.nMax)," with gcd(N,6E)=1. Neither bound is asserted for the original problem."));
    PPLine(stream,"The N rows display the order test L divides N and the product test D_E divides uN, even when they fail. No M, h or p is chosen. Conditions requiring those choices are not tested.");
    PPLine(stream,"Here L=ord_a(q), computed as the lcm of the prime-power orders. A nonunit gcd(L,6E) excludes every N with gcd(N,6E)=1 for that row.");
    PPLine(stream,"S_ℓ consists of primes s dividing j=(ℓ−1)/ord_ℓ(2) with s not dividing 6E. For every N with gcd(N,6E)=1, the full factor ℓ^f is omitted from D_E exactly when some s in S_ℓ divides N. Empty S_ℓ means that ℓ^f is always retained. This description has no upper bound on N.");
    PPLine(stream,"The Frobenius tables include all s occurring in the N evaluations and every s in every S_ℓ, including those above the N range.");
    PPLine(stream,"An ordinary odd prime divisor ℓ of q−1 means v_ℓ(2^ord_ℓ(2)−1)=1. Such a prime cannot divide N in a solution. The column records this test separately.");
    PPLine(stream,"The indices d in each factor table distinguish the primitive exponent t=ord_ℓ(2) from other cyclotomic values Φ_d(2) containing ℓ. The q factor itself is removed completely; its primitive exponent is r.");
    PPLine(stream,"Factorization products and the primality of every factor are checked. Above 2^64, a GAP primality proof is required.");
    PPLine(stream,"\n## Frobenius and the selected factors\n");
    PPLine(stream,"For ℓ not dividing m, arithmetic Frobenius in Q(ζ_m) sends ζ_m to ζ_m^ℓ; its order is ord_m(ℓ), the residue degree. See [Milne, Algebraic Number Theory, cyclotomic extensions and Frobenius](https://www.jmilne.org/math/CourseNotes/ANT.pdf).");
    PPLine(stream,"For t=ord_ℓ(2), ℓ=1 modulo t, so Frobenius in Q(ζ_t) is the identity. To form D_E, use Q(ζ_s) for each prime s dividing N. Since gcd(E,N)=1 and t divides E, gcd((ℓ−1)/t,N)=gcd(ℓ−1,N).");
    PPLine(stream,"Thus ℓ^f belongs to D_E precisely when Frobenius is nonidentity in every Q(ζ_s) with s dividing N and s≠ℓ. At s=ℓ the prime is ramified; it imposes no exclusion from D_E, since ℓ does not divide ℓ−1. The table retains the full ℓ^f, and checks D_E divides ((q−1)/2^c)N.");
    PPLine(stream,"For ℓ≠s, Φ_s over F_ℓ has (s−1)/f distinct irreducible factors of degree f=ord_s(ℓ); these correspond to the prime ideals above ℓ. At ℓ=s, Φ_s=(X−1)^(s−1) over F_s. The displayed decomposition data can therefore be checked by finite-field factorization.");
    PPLine(stream,"The order in Q(ζ_N) alone is insufficient to select D_E: the restrictions to each Q(ζ_s) matter. This is the exact role of the additional cyclotomic moduli here.");
    PPLine(stream,"\n## Reproduction\n\n```sh\n~/fluid/sauce/gap/gap -A -q -b mersenne_tables.g\n```\n");
    PPLine(stream,"To change the bounds, in GAP:\n\n```gap\nPPConfig := rec(qMax := 10^8, nMax := 121, factorE := 256,\n                output := \"mersenne-tables.md\");;\nRead(\"mersenne_tables.g\");\n```\n");
    PPLine(stream,"`mersenne_tables_test.g` checks the arithmetic and table helpers on independent small fixtures; those fixtures are not included as candidates.");
    CloseStream(stream);
    Print("Wrote ",config.output,"; counts = ",scan.counts,"\n");
    return scan;
end;

if not IsBound(PPLibraryOnly) or not PPLibraryOnly then
    if not IsBound(PPConfig) then
        PPConfig := rec(qMax := 10^8,nMax := 121,factorE := 256,
            output := "mersenne-tables.md");
    fi;
    PPResult := PPReport(PPConfig);;
fi;
