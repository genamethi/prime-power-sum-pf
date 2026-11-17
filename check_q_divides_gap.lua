-- Finite proof check for the branch q | M in prime_power_notes.tex.
-- Analytic inputs proved there: q < 10^6 and M < 4 * 10^6.
-- Every rejection applies to all M, N and p satisfying the gap equations.
-- Run with Lua 5.5, Lua 5.4, or LuaJIT; no external libraries are needed.

local Q_LIMIT = 1000000
local M_LIMIT = 4000000

local function filter(values, predicate)
    local result = {}
    for _, value in ipairs(values) do
        if predicate(value) then result[#result + 1] = value end
    end
    return result
end

local function primes_below(limit)
    local composite, primes = {}, {}
    for p = 2, math.floor(math.sqrt(limit - 1)) do
        if not composite[p] then
            for n = p * p, limit - 1, p do composite[n] = true end
        end
    end
    for n = 2, limit - 1 do
        if not composite[n] then primes[#primes + 1] = n end
    end
    return primes
end

local function valuation(n, prime)
    local exponent = 0
    while n % prime == 0 do
        n = math.floor(n / prime)
        exponent = exponent + 1
    end
    return exponent
end

local function gcd(a, b)
    while b ~= 0 do a, b = b, a % b end
    return a
end

local function lcm(a, b)
    return (a / gcd(a, b)) * b
end

local function power_modulo(base, exponent, modulus)
    local result = 1
    base = base % modulus
    while exponent > 0 do
        if exponent % 2 == 1 then result = result * base % modulus end
        base = base * base % modulus
        exponent = math.floor(exponent / 2)
    end
    return result
end

local function prime_factors(n, primes)
    local factors = {}
    for _, prime in ipairs(primes) do
        if prime * prime > n then break end
        if n % prime == 0 then
            factors[#factors + 1] = prime
            repeat n = math.floor(n / prime) until n % prime ~= 0
        end
    end
    if n > 1 then factors[#factors + 1] = n end
    return factors
end

local function order_of_two(prime, primes)
    local order = prime - 1
    for _, factor in ipairs(prime_factors(order, primes)) do
        while order % factor == 0
            and power_modulo(2, order / factor, prime) == 1 do
            order = math.floor(order / factor)
        end
    end
    return order
end

local function meets_base_restrictions(q)
    return q > 3
        and q % 3 == 1
        and valuation(q - 1, 2) % 2 == 0
        and power_modulo(2, (q - 1) / 3, q) == 1
end

local function forced_gap_divisor(q, primes)
    -- r | M, q | M, 2 | M, and v3(M) = v3(q-1) - 1.
    local r = order_of_two(q, primes)
    local three_part = 3 ^ (valuation(q - 1, 3) - 1)
    local L = q * lcm(lcm(2, r), three_part)
    if L >= M_LIMIT then return L end

    -- The odd part of q-1 divides A. Each of its primes therefore
    -- forces its order at base 2 to divide M as well.
    for _, prime in ipairs(prime_factors(q - 1, primes)) do
        if prime > 2 then
            L = lcm(L, order_of_two(prime, primes))
            if L >= M_LIMIT then return L end
        end
    end
    return L
end

local witnesses = {
    -- ell, ord_ell(2), largest divisor of ell-1 coprime to 6.
    { prime = 5,  two_order = 4, allowed_q_order = 1 },
    { prime = 7,  two_order = 3, allowed_q_order = 1 },
    { prime = 31, two_order = 5, allowed_q_order = 5 },
}

local function first_obstruction(q, primes)
    local L = forced_gap_divisor(q, primes)
    if L >= M_LIMIT then return "minimum M exceeds height bound" end
    if valuation(L, 3) > valuation(q - 1, 3) - 1 then
        return "incompatible valuation at 3"
    end

    -- If ord_ell(2) | L and ell ~= q, then ell | A.
    -- Since gcd(N,6)=1, the order of q modulo ell must divide
    -- the largest divisor of ell-1 coprime to 6.
    for _, witness in ipairs(witnesses) do
        if q ~= witness.prime and L % witness.two_order == 0
            and power_modulo(q, witness.allowed_q_order, witness.prime) ~= 1 then
            return "incompatible order modulo " .. witness.prime
        end
    end
    return nil
end

local function main()
    local primes = primes_below(Q_LIMIT)
    local candidates = filter(primes, meets_base_restrictions)
    local survivors, counts = {}, {}
    for _, q in ipairs(candidates) do
        local obstruction = first_obstruction(q, primes)
        if obstruction then
            counts[obstruction] = (counts[obstruction] or 0) + 1
        else
            survivors[#survivors + 1] = q
        end
    end

    print("Primes satisfying the base restrictions: " .. #candidates)
    for _, reason in ipairs({
        "minimum M exceeds height bound",
        "incompatible valuation at 3",
        "incompatible order modulo 5",
        "incompatible order modulo 7",
        "incompatible order modulo 31",
    }) do
        print(reason .. ": " .. (counts[reason] or 0))
    end
    print("Survivors: " .. #survivors)
    for _, q in ipairs(survivors) do print(q) end
    assert(#survivors == 0, "The stated restrictions did not exclude every q")
end

-- Exact arithmetic also on LuaJIT: modular products are below 10^12.
-- The initial L is below 10^12; subsequent lcm products are below
-- 4 * 10^12 because we stop as soon as L reaches M_LIMIT.
-- All integer intermediates are therefore below 2^53.
main()
