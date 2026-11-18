-- Exact finite check for M = 2^a q^b, using the external bound q < 10^6.
-- Run with: lua check_power_two_orders.lua

local LIMIT = 1000000

local function filter(values, predicate)
    local result = {}
    for _, value in ipairs(values) do
        if predicate(value) then
            result[#result + 1] = value
        end
    end
    return result
end

local function primes_below(limit)
    -- Sieve of Eratosthenes; mutation stays local to the sieve.
    local composite = {}
    for p = 2, math.floor(math.sqrt(limit - 1)) do
        if not composite[p] then
            for multiple = p * p, limit - 1, p do
                composite[multiple] = true
            end
        end
    end

    local primes = {}
    for n = 2, limit - 1 do
        if not composite[n] then
            primes[#primes + 1] = n
        end
    end
    return primes
end

local function valuation_at_two(n)
    local exponent = 0
    while n % 2 == 0 do
        n = math.floor(n / 2)
        exponent = exponent + 1
    end
    return exponent
end

local function power_modulo(base, exponent, modulus)
    -- Repeated squaring. With modulus < 10^6, products stay below 10^12,
    -- so the arithmetic is exact even in Lua versions using doubles.
    local result = 1
    base = base % modulus
    while exponent > 0 do
        if exponent % 2 == 1 then
            result = (result * base) % modulus
        end
        base = (base * base) % modulus
        exponent = math.floor(exponent / 2)
    end
    return result
end

local function is_one_modulo_three(q)
    return q % 3 == 1
end

local function has_even_two_valuation(q)
    return valuation_at_two(q - 1) % 2 == 0
end

local function has_power_of_two_order(q)
    -- For prime q, ord_q(2) divides q - 1. It is a power of two
    -- exactly when it divides the largest power of two dividing q - 1.
    local two_part = 2 ^ valuation_at_two(q - 1)
    return power_modulo(2, two_part, q) == 1
end

local function is_one_modulo_five(q)
    return q % 5 == 1
end

local function display(label, values)
    local entries = {}
    for i, value in ipairs(values) do
        entries[i] = string.format("%.0f", value)
    end
    print(label .. ": [" .. table.concat(entries, ", ") .. "]")
end

local function main()
    -- Apply the restrictions in the order of the mathematical argument.
    local primes = primes_below(LIMIT)
    local one_modulo_three = filter(primes, is_one_modulo_three)
    local even_valuation = filter(one_modulo_three, has_even_two_valuation)
    local candidates = filter(even_valuation, has_power_of_two_order)
    local survivors = filter(candidates, is_one_modulo_five)

    display("Power-of-two order, q = 1 mod 3, even v2(q-1)", candidates)
    display("Also q = 1 mod 5", survivors)

    assert(#candidates == 1 and candidates[1] == 274177,
        "Expected exactly the candidate 274177")
    assert(#survivors == 0, "Expected no survivors")
end

main()
