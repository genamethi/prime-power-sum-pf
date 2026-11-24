\\ From the project root: gp -q scripts/check_odd_part.gp
\\ Each row gives a prime and its exponent in 2^(256*b)-1.
odd_part_factors(b) = {
  my(primes = [3, 5, 17, 257, 65537, 641, 6700417, 274177,
               67280421310721, 59649589127497217, 5704689200685129054721]);
  my(extra = [7, 3, 13, 241, 97, 673, 193, 22253377,
              18446744069414584321, 769, 442499826945303593556473164314770689,
              349621839326921795694385454593,
              331192380488114152600457428497953408512758882817]);
  if (b == 3, primes = concat(primes, extra), if (b != 1, error("use b=1 or 3")));
  my(factors = factorback_checked(primes, b));
  return(factors);
}

factorback_checked(primes, b) = {
  my(distinct = Set(primes), k = #distinct);
  my(factors = matrix(k, 2));
  for (i = 1, k,
    my(p = distinct[i]);
    if (!isprime(p), error("composite entry: ", p));
    factors[i, 1] = p;
    factors[i, 2] = sum(j = 1, #primes, primes[j] == p);
    if (p > 3 && factors[i, 2] > 1, error("repeated possible q: ", p));
  );
  if (factorback(factors) != 2^(256*b)-1, error("wrong product for ", b));
  return(factors);
}

possible_u(factors, N, modulus) = {
  my(residues = [1]);
  for (i = 1, matsize(factors)[1],
    my(p = factors[i, 1], exponent = factors[i, 2] - valuation(N, p));
    my(including = []);
    if (exponent > 0,
      including = Set(apply(u -> (u*p^exponent) % modulus, residues)));
    if (gcd(N, p-1) == 1,
      residues = including,
      residues = setunion(residues, including));
    if (#residues == 0, break);
  );
  return(residues);
}

remaining_gaps(factors, bits) = {
  my(modulus = 2^bits, gaps = List());
  forstep (N = 5, 2519, 2,
    if (gcd(N, 30) != 1, next);
    my(required = lift(-1/Mod(N, modulus)));
    if (setsearch(possible_u(factors, N, modulus), required), listput(gaps, N));
  );
  return(Vec(gaps));
}

{
  my(shared = List());
  forprime (p = 3, 2519,
    if (Mod(2, p^2)^(p-1) == 1, listput(shared, p)));
  if (Vec(shared) != [1093], error("unexpected valuation list"));
  print("possible common primes of N and q-1: ", Vec(shared));
  forstep (b = 1, 3, 2,
    my(factors = odd_part_factors(b));
    for (bits = 8, 10,
      my(gaps = remaining_gaps(factors, bits));
      print("odd part ", b, ", modulus ", 2^bits, ": ", gaps);
      if (bits == 10 && #gaps, error("exclusion failed")));
  );
}
