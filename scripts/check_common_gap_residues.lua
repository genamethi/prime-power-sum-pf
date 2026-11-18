-- Necessary residue conditions for d | gcd(M,N), with d prime.
-- Inputs from prime_power_notes.tex: p = 1 mod d, n1 even,
-- M = hE with E even, d | h, and (2h-1)n1 < N.
-- The modulo-5 split retains the exponent classes modulo 4.
-- This checks whole residue classes, not primes q below a search bound.
-- Lua 5.5, 5.4 and LuaJIT; no external libraries.

local function gcd(a, b)
    while b ~= 0 do a, b = b, a % b end
    return a
end

local function lcm(a, b)
    return math.floor(a / gcd(a, b)) * b
end

local function power_modulo(base, exponent, modulus)
    local result = 1
    while exponent > 0 do
        if exponent % 2 == 1 then result = result * base % modulus end
        base = base * base % modulus
        exponent = math.floor(exponent / 2)
    end
    return result
end

local function primes_through(limit)
    local result = {}
    for n = 2, limit do
        local prime = true
        for divisor = 2, math.floor(math.sqrt(n)) do
            if n % divisor == 0 then prime = false; break end
        end
        if prime then result[#result + 1] = n end
    end
    return result
end

local function even_exponents_by_value(d, period)
    local result = {}
    for exponent = 0, period - 1, 2 do
        local value = power_modulo(2, exponent, d)
        result[value] = result[value] or {}
        result[value][#result[value] + 1] = exponent
    end
    return result
end

local function passes_modulo_five(c, m1, e, N)
    if (m1 - c) % 4 == 0 then
        return c % 4 == 0 and N % 5 ~= 0
    end
    return (c + e) % 4 == 2
end

local function possible_even_E(c, m1, h, period)
    -- hE = m1-c mod period is soluble with E even exactly in this case.
    return (m1 - c) % gcd(2 * h, period) == 0
end

local function local_residues(d, N, e, h)
    local period = lcm(4, d - 1)
    local exponents = even_exponents_by_value(d, period)
    local result = {}
    for q = 2, d - 1 do
        -- q = 0,1 mod d are already excluded by gcd(g,q(q-1))=1.
        local small_power = (1 - power_modulo(q, N + e, d)) % d
        local large_power = (1 - power_modulo(q, e, d)) % d
        local found = false
        for _, c in ipairs(exponents[small_power] or {}) do
            for _, m1 in ipairs(exponents[large_power] or {}) do
                if passes_modulo_five(c, m1, e, N)
                    and possible_even_E(c, m1, h, period) then
                    found = true
                    break
                end
            end
            if found then break end
        end
        if found then result[#result + 1] = q end
    end
    return result
end

local function surviving_families(d, N)
    local result = {}
    local e = 2
    while (2 * d - 1) * e < N do
        local h = d
        while (2 * h - 1) * e < N do
            local residues = local_residues(d, N, e, h)
            if #residues > 0 then
                result[#result + 1] = {e = e, h = h, residues = residues}
            end
            h = h + d
        end
        e = e + 2
    end
    return result
end

local function check_global_seventeen()
    -- All relevant exponents are periodic modulo 16 (also a multiple of 4).
    -- Allow the 4 | M branch regardless of 5 | N. Ignoring 3 | N
    -- and allowing h=17 also enlarge the necessary residue set.
    for e = 2, 16, 2 do
        for N = 1, 15, 2 do
            local representative = N % 5 == 0 and N + 16 or N
            assert(#local_residues(17, representative, e, 17) == 0,
                "Unexpected residue class for 17 | gcd(M,N)")
        end
    end
end

local function main()
    check_global_seventeen()
    print("17 | gcd(M,N): no residue classes")

    local excluded = {}
    for _, d in ipairs(primes_through(200)) do
        if d >= 7 and #surviving_families(d, 5 * d) == 0 then
            excluded[#excluded + 1] = d
        end
    end
    print("N=5d forces coprime gaps for d in: " .. table.concat(excluded, ", "))

    for _, pair in ipairs({{49, 7}, {55, 11}, {65, 13}, {77, 7},
        {77, 11}, {85, 17}, {91, 7}, {91, 13}, {121, 11}}) do
        local N, d = pair[1], pair[2]
        local families = surviving_families(d, N)
        print(string.format("N=%d, d=%d: %d local families", N, d, #families))
        for _, family in ipairs(families) do
            print(string.format("  n1=%d, h=%d, q mod d in {%s}",
                family.e, family.h, table.concat(family.residues, ",")))
        end
    end
end

main()
