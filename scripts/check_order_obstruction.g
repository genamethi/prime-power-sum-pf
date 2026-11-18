# Exact certificates excluding all N with gcd(N,E)=1 for a fixed q.
# Run: ~/fluid/sauce/gap/gap -A -q -b --quitonbreak check_order_obstruction.g
# This wider experiment does NOT impose the remaining-case conditions on q.
PPLibraryOnly := true;;
Read("mersenne_tables.g");
Read("mersenne-factors-256.g");

PPOrderFromMultiple := function(base,prime,multiple)
    local order,s;
    order:=multiple;
    if PowerModInt(base,order,prime)<>1 then Error("not an order multiple"); fi;
    for s in Set(FactorsInt(multiple)) do
        while order mod s=0 and PowerModInt(base,QuoInt(order,s),prime)=1 do
            order:=QuoInt(order,s);
        od;
    od;
    return order;
end;

PPOrderWitness := function(q,r,factors)
    local pair,ell,f,modulus,phi,s,v,w,residue;
    for pair in factors do
        ell := pair[1]; f := pair[2];
        if ell=q then continue; fi;
        modulus := ell^f;
        phi := (ell-1)*ell^(f-1);
        for s in Set(FactorsInt(r)) do
            v := PPValuation(phi,s);
            if v=0 then continue; fi;
            w := QuoInt(phi,s^v);
            residue := PowerModInt(q,w,modulus);
            if residue<>1 then
                return rec(ell:=ell,f:=f,s:=s,w:=w,residue:=residue);
            fi;
        od;
    od;
    return fail;
end;

PPOrderExperiment := function(limit,output)
    local stream,E,factors,pair,q,r,rho,witness,count,excluded,missing,evenRho;
    stream := OutputTextFile(output,false);
    SetPrintFormattingStatus(stream,false);
    PPLine(stream,"# An order obstruction for all N\n");
    PPLine(stream,"Put r=ord_q(2), E=lcm(2,r), rho=v_q(2^r−1), a=(2^E−1)/q^rho. In the remaining case, a divides A, ord_a(q) divides N, and gcd(E,N)=1.");
    PPLine(stream,"\n## Proven certificate\n");
    PPLine(stream,"Let ℓ^f divide a, let s be a prime divisor of r, and write φ(ℓ^f)=s^v w with s not dividing w. If q^w is not 1 modulo ℓ^f, then s divides ord_(ℓ^f)(q). Consequently s divides both E and N, a contradiction. This excludes every N for that q.");
    PPLine(stream,"Proof: ord_(ℓ^f)(q) divides s^v w. If its order were coprime to s, it would divide w, contrary to the displayed residue.");
    PPLine(stream,"\n## Research question\n");
    PPLine(stream,"Can the factorization of 2^E−1 force such a certificate? Equivalently, can one prove gcd(r,ord_a(q))>1, either uniformly under our remaining assumptions or outside an effectively bounded set? This is not proved here.");
    PPLine(stream,"The experiment below first tests the stronger question for every odd prime q with E in the stated range, without imposing the q residue conditions or even rho. It distinguishes evidence for that wider question from the repeated-factor case required by the notes.");
    PPLine(stream,"\n## Kummer formulation\n");
    PPLine(stream,"For a prime-level certificate, let s divide r and ℓ−1, v=v_s(ℓ−1), w=(ℓ−1)/s^v, K=Q(ζ_(s^v)), and F=K(q^(1/s^v)). A prime P of K above ℓ has residue field F_ℓ. For α=q^(1/s^v), arithmetic Frobenius satisfies Frob_P(α)/α = q^w modulo a prime of F above P. Hence q^w≠1 detects nontrivial Frobenius in F/K and forces s to divide ord_ℓ(q). See [Stevenhagen, Kummer Theory and Reciprocity Laws](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-Kummer.pdf), Sections 2 and 4.");
    PPLine(stream,"The task is to force this nontrivial Frobenius at a prime dividing the specific integer a. If the witness only appears modulo ℓ^f with f>1 and s=ℓ, use the prime-power order test above; the prime-level Kummer formulation does not cover that case.");
    PPLine(stream,"\n## Exact finite experiment\n");
    PPLine(stream,Concatenation("Enumerate every even E with 2 ≤ E ≤ ",String(limit),
        ", factor 2^E−1, and retain precisely the prime factors q with lcm(2,ord_q(2))=E. Each q therefore occurs once. The supplied factorizations in `mersenne-factors-256.g` were computed with PARI/GP; GAP verifies their products and proves primality of every factor. The residue column is an exact modular-power certificate; no bound on N or p is used.\n"));
    PPRow(stream,["q","r","E","rho","ℓ","f","s dividing r","w","q^w mod ℓ^f"]);
    PPRow(stream,List([1..9],x -> "---"));
    count:=0; excluded:=0; missing:=[]; evenRho:=0;
    for E in [2,4..limit] do
        factors := PPMersenneFactorCache[E];
        if Product(factors,pair -> pair[1]^pair[2])<>2^E-1 or
            not ForAll(factors,pair -> PPProvenPrime(pair[1])) then
            Error("uncertified cached factorization at E=",E);
        fi;
        for pair in factors do
            q:=pair[1]; rho:=pair[2]; r:=PPOrderFromMultiple(2,q,E);
            if Lcm(2,r)<>E then continue; fi;
            count:=count+1;
            if rho mod 2=0 then evenRho:=evenRho+1; fi;
            witness:=PPOrderWitness(q,r,factors);
            if witness=fail then
                Add(missing,[q,r,E,rho,QuoInt(2^E-1,q^rho)]);
            else
                excluded:=excluded+1;
                PPRow(stream,[q,r,E,rho,witness.ell,witness.f,witness.s,witness.w,witness.residue]);
            fi;
        od;
    od;
    PPLine(stream,"\n## Results\n");
    PPLine(stream,Concatenation(String(count)," prime choices; ",String(excluded),
        " certified order obstructions. Cases without a certificate, as [q,r,E,rho,a]: ",String(missing),"."));
    PPLine(stream,Concatenation("Number with even rho: ",String(evenRho),". The experiment does not establish the proposed uniform statement."));
    PPLine(stream,"\nReproduce with:\n\n```sh\n~/fluid/sauce/gap/gap -A -q -b --quitonbreak check_order_obstruction.g\n```\n");
    CloseStream(stream);
    Print(count," choices; ",excluded," witnesses; missing = ",missing,
        "; even rho = ",evenRho,"\n");
end;

PPOrderExperiment(256,"order-obstructions.md");
