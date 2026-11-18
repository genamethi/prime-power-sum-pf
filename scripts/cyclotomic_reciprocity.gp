\\ Relative norms and power-residue symbols; PARI/GP, no optional packages.
\\ Load with read("cyclotomic_reciprocity.gp").
\\ Examples: reciprocity_case(1777,37,37), reciprocity_case(1093,364,13).

require(condition, message) = if(!condition, error(message));

power_valuation(base, exponent, prime) = {
  my(v = 0);
  while(Mod(base, prime^(v+1))^exponent == 1, v++);
  v;
};

\\ Norm from Q(zeta_t) to Q(zeta_m), with zeta_m=zeta_t^(t/m).
\\ Inclusion-exclusion on the exponents coprime to t.
relative_cyclotomic_norm(t, m) = {
  require(t % m == 0 && m > 1, "m must divide t");
  my(z = Mod(x, polcyclo(m)), h = t/m, result = 1);
  fordiv(h, d,
    if(gcd(d,m) == 1 && moebius(d),
      result *= (2^(h/d) - z^lift(Mod(d,m)^(-1)))^moebius(d)
    )
  );
  result = lift(result);
  require(denominator(result) == 1, "norm must be integral");
  require(polresultant(polcyclo(m),result) == subst(polcyclo(t),x,2),
          "absolute norm mismatch");
  result;
};

\\ The logarithm is only in the small group of m-th roots of unity.
root_exponent(value, root, m) = {
  my(power = root^0);
  for(j=0,m-1,
    if(value == power, return(j));
    power *= root
  );
  error("value is not an m-th root of unity");
};

\\ Product of the tame local symbols (pi,q)_m at all primes above q.
\\ Rows: [embedding i, v_P(pi), unit modulo q, exponent relative to zeta_m].
polynomial_symbols_above_q(q, pi, m) = {
  require(isprime(q) && (q-1)%m == 0, "q must split in Q(zeta_m)");
  my(rho = valuation(polresultant(polcyclo(m),pi),q), precision = rho+2);
  my(z = znprimroot(q)^((q-1)/m));
  my(zlift = Mod(lift(z),q^precision)^(q^(precision-1)));
  my(rows = List(), total = 0);
  require(zlift^m == 1, "root lift failed");
  for(i=1,m-1, if(gcd(i,m)==1,
    my(value = lift(subst(pi,x,zlift^i)));
    require(value != 0, "insufficient q-adic precision");
    my(v = valuation(value,q), unit = Mod(value/q^v,q));
    my(exponent = root_exponent(((-1)^v*unit)^((q-1)/m),z^i,m));
    listput(rows,[i,v,lift(unit),exponent]);
    total = (total+exponent)%m
  ));
  [pi,Vec(rows),total];
};

symbols_above_q(q,t,m) = polynomial_symbols_above_q(q,relative_cyclotomic_norm(t,m),m);

\\ Independent ideal factorization and residue-field computation.
\\ Omit primes above q and above the prime dividing m.
\\ Rows: [ell, residue degree, ideal valuation, symbol exponent].
symbols_away_from_q(q,t,m) = {
  my(pi = relative_cyclotomic_norm(t,m), nf = nfinit(polcyclo(m)));
  my(factors = idealfactor(nf,pi), rows = List(), total = 0);
  for(j=1,matsize(factors)[1],
    my(P=factors[j,1], multiplicity=factors[j,2], ell=P.p);
    if(ell == q || m%ell == 0, next);
    my(reduction = nfmodprinit(nf,P));
    my(z = nfmodpr(nf,x,reduction), Q = ell^P.f);
    my(value = nfmodpr(nf,q,reduction)^((Q-1)/m));
    my(exponent = root_exponent(value,z,m));
    listput(rows,[ell,P.f,multiplicity,exponent]);
    total = (total + multiplicity*exponent)%m
  );
  [Vec(rows),total];
};

\\ Sufficient condition for the ramified local symbols to vanish.
\\ At 2, q congruent to 1 mod 2m gives an unramified Kummer extension;
\\ pi is a 2-adic unit. For m=2, q>0 also handles the real place.
local_correction_vanishes(q,m) = {
  my(f=factor(m));
  require(matsize(f)[1] == 1, "m must be a prime power");
  my(s=f[1,1], b=f[1,2]);
  if(s==2, return(q%(2*m)==1));
  q%(s*m)==1;
};

reciprocity_case(q,t,m) = {
  my(atq=symbols_above_q(q,t,m), away=symbols_away_from_q(q,t,m));
  my(vanishes=local_correction_vanishes(q,m));
  if(vanishes, require(atq[3]==away[2], "reciprocity check failed"));
  print("q=",q,", t=",t,", m=",m,", pi=",atq[1]);
  print("above q: ",atq[2]);
  print("away: ",away[1]);
  print("exponents: above q=",atq[3],", away=",away[2],
        "; local correction certified trivial=",vanishes);
  [atq,away,vanishes];
};

\\ Group all cyclotomic factors whose indices have the full s-part of E.
\\ The absolute norm is H=(2^E-1)/(2^(E/s)-1), for every m=s^b dividing r.
\\ Return a residue modulo q: nontriviality is independent of the root chosen.
grouped_residue(q,m) = {
  my(r=znorder(Mod(2,q)), f=factor(m));
  require(matsize(f)[1]==1 && r%m==0, "m must be a prime power dividing r");
  my(s=f[1,1], E=lcm(2,r), B=2^(E/m), H=(2^E-1)/(2^(E/s)-1));
  my(rho=valuation(H,q), cofactor=H/q^rho, z=Mod(B,q));
  my(derivative=subst(deriv(polcyclo(m)),x,z));
  my(W=(-1)^rho*Mod(cofactor,q)/derivative);
  for(i=2,m-1, if(gcd(i,m)==1,
    W *= (z-z^i)^lift(Mod(i,m)^(-1))
  ));
  my(residue=W^((q-1)/m), exponent=root_exponent(residue,z,m));
  [r,rho,B,H,cofactor,lift(W),lift(residue),exponent];
};

\\ An independent computation from a complete, certified factorization of H.
\\ Each unramified factor gives one prime ideal, zeta_m reducing to B.
grouped_factor_symbols(q,m,factors) = {
  my(data=grouped_residue(q,m), B=data[3], H=data[4]);
  require(prod(j=1,matsize(factors)[1],factors[j,1]^factors[j,2])==H,
          "factorization does not equal H");
  my(rows=List(), total=0);
  for(j=1,matsize(factors)[1],
    my(ell=factors[j,1], multiplicity=factors[j,2]);
    require(isprime(ell), "uncertified prime factor");
    if(ell==q || m%ell==0, next);
    require((ell-1)%m==0 && znorder(Mod(B,ell))==m, "wrong residue degree");
    my(residue=Mod(q,ell)^((ell-1)/m));
    my(exponent=root_exponent(residue,Mod(B,ell),m));
    listput(rows,[ell,multiplicity,lift(residue),exponent]);
    total=(total+multiplicity*exponent)%m
  );
  [Vec(rows),total];
};

check_grouped_case(q,m,{factors=[]}) = {
  my(data=grouped_residue(q,m), atq=polynomial_symbols_above_q(q,data[3]-x,m));
  require(atq[3]==data[8], "q-adic and derivative tests disagree");
  my(away=[]);
  if(#factors,
    away=grouped_factor_symbols(q,m,factors);
    if(local_correction_vanishes(q,m),
      require(away[2]==data[8], "grouped reciprocity failed")
    )
  );
  [data,atq[2],away,local_correction_vanishes(q,m)];
};
