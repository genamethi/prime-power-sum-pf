# Run: ~/fluid/sauce/gap/gap -A -q -b --quitonbreak mersenne_tables_test.g
PPLibraryOnly := true;;
Read("mersenne_tables.g");

PPTestCount := 0;;
PPCheck := function(condition,message)
    if not condition then Error(message); fi;
    PPTestCount := PPTestCount+1;
end;

# Compare modular valuations with actual integer powers.
for prime in [3,5,7,11,13,1093,3511] do
    r := OrderMod(2,prime);
    PPCheck(PPPowerValuation(2,r,prime)=PPValuation(2^r-1,prime),
        "valuation mismatch");
od;
PPCheck(PPPowerValuation(2,364,1093)=2,"Wieferich valuation fixture");
PPCheck(not PPProvenPrime(341),"composite fixture");
PPCheck(PPProvenPrime(2^89-1),"large certified prime fixture");

# Include all three residue ranges, small limits, and inclusive endpoints.
PPCheck(IsEmpty(PPStructuralScan(rec(qMax:=96,factorE:=256)).rows),"lower endpoint");
PPCheck(IsEmpty(PPStructuralScan(rec(qMax:=432,factorE:=256)).rows),"upper endpoint");
PPCheck(List(PPStructuralScan(rec(qMax:=433,factorE:=256)).rows,b -> b.q)=[433],
    "inclusive prime endpoint");
scan := PPStructuralScan(rec(qMax:=10000,factorE:=256));;
PPCheck(List(scan.rows,b -> b.q)=[433,1777,2113,2833,4177],"bounded structural scan");

# A deliberately nonadmissible fixture checks that no row is relabelled a solution.
PPCheck(not PPBaseCompatible(PPBaseData(433)),"odd rho must fail the full tests");
PPCheck(IsEmpty(PPScan(10000).candidates),"strict scan must remain empty");

for b in scan.rows do
    fs := PPFactorsOfA(b);
    a := QuoInt(2^b.E-1,b.q^b.e);
    PPCheck(Product(fs,f -> f.ell^f.f)=a,"factor product");
    PPCheck(Gcd(a,b.q)=1,"q removed completely");
    for f in fs do
        modulus := f.ell^f.f;
        PPCheck(PowerModInt(b.q,f.orderQPower,modulus)=1,"prime-power order");
        for divisor in Filtered(Set(FactorsInt(f.orderQPower)),p -> p>1) do
            PPCheck(PowerModInt(b.q,QuoInt(f.orderQPower,divisor),modulus)<>1,
                "order minimality");
        od;
        for d in f.phiIndices do
            quotient := QuoInt(d,f.t);
            PPCheck(d mod f.t=0 and quotient=f.ell^PPValuation(quotient,f.ell),
                "cyclotomic occurrence must have index t*ell^j");
        od;
    od;
    for N in Filtered([5..121],N -> Gcd(N,6*b.E)=1) do
        # Independently select factors by Frobenius on each prime cyclotomic field.
        selected := Filtered(fs,f -> ForAll(Set(FactorsInt(N)),s -> f.ell mod s<>1));
        PPCheck(PPSelectedProduct(fs,N)=Product(selected,f -> f.ell^f.f),
            "Frobenius selection and residual-index selection disagree");
    od;
od;

# Reduce the integer cyclotomic polynomial, then factor over the finite field.
# GAP's CyclotomicPolynomial(GF(s),s) removes inseparability: it is not the
# reduction of the integer polynomial needed for the ramified check.
PPCheckPolynomialFactorizations := function()
    local ell,s,factors,frob,coefficients;
    for ell in [3,5,7,11,13,17,19,31] do
      for s in [5,7,11] do
        coefficients := CoefficientsOfUnivariatePolynomial(CyclotomicPolynomial(Rationals,s));
        factors := Factors(UnivariatePolynomial(GF(ell),
            List(coefficients,c -> c*One(GF(ell)))));
        frob := PPFrobenius(ell,s);
        if ell=s then
            PPCheck(frob[1]="ramified" and Length(factors)=s-1 and
                Length(Set(factors))=1 and Degree(factors[1])=1,"ramification");
        else
            PPCheck(List(factors,Degree)=List([1..QuoInt(s-1,frob[2])],i -> frob[2]),
                "residue degrees versus finite-field factorization");
            PPCheck(GaloisCyc(E(s),ell)=E(s)^frob[1],"cyclotomic automorphism");
        fi;
      od;
    od;
end;
PPCheckPolynomialFactorizations();
Print("Passed ",PPTestCount," checks.\n");
QUIT;
