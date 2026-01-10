"""
Common benchmark definitions shared between FFTA and FFTW benchmarks
"""

using Primes

# Cumulative products of 3,4,5,5,7,11
const CUMULATIVE_PRODUCTS = cumprod([3, 4, 5, 5, 7, 11])

# All primes below 20000
const PRIMES_BELOW_20000 = primes(20000)

# Odd powers of 2
const ODD_POWERS_OF_2 = [2^i for i in 1:2:15]

# Even powers of 2
const EVEN_POWERS_OF_2 = [2^i for i in 2:2:14]

# Powers of 3
const POWERS_OF_3 = [3^i for i in 1:9]

"""
    create_size_categories()

Create a dictionary mapping array sizes to their categories.
"""
function create_size_categories()
    categories = Dict{Int, String}()

    # Categorize odd powers of 2
    for n in ODD_POWERS_OF_2
        categories[n] = "odd_power_of_2"
    end

    # Categorize even powers of 2
    for n in EVEN_POWERS_OF_2
        categories[n] = "even_power_of_2"
    end

    # Categorize powers of 3
    for n in POWERS_OF_3
        categories[n] = "power_of_3"
    end

    # Categorize cumulative products
    for n in CUMULATIVE_PRODUCTS
        categories[n] = "cumulative_product"
    end

    # Categorize primes
    for n in PRIMES_BELOW_20000
        categories[n] = "prime"
    end

    return categories
end

"""
    get_all_sizes()

Get all benchmark sizes in a flat array.
"""
function get_all_sizes()
    return vcat(ODD_POWERS_OF_2, EVEN_POWERS_OF_2, POWERS_OF_3,
                CUMULATIVE_PRODUCTS, PRIMES_BELOW_20000)
end
