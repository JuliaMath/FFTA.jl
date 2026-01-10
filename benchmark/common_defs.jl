"""
Common benchmark definitions shared between FFTA and FFTW benchmarks
"""

using Primes

# Cumulative products of 3,4,5,5,7,11
const CUMULATIVE_PRODUCTS = cumprod([3, 4, 5, 5, 7, 11])

# Select 20 primes with logarithmic spacing
const SELECTED_PRIMES = begin
    # Generate logarithmically spaced target points
    all_primes = primes(20000)
    log_points = exp.(range(log(2), log(20000), length=20))

    # Find nearest prime to each logarithmic point
    selected = Int[]
    for target in log_points
        # Find the prime closest to this target
        idx = argmin(abs.(all_primes .- target))
        candidate = all_primes[idx]
        # Avoid duplicates
        if !(candidate in selected)
            push!(selected, candidate)
        end
    end

    sort!(selected)
end

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

    # Categorize composite numbers
    for n in CUMULATIVE_PRODUCTS
        categories[n] = "composite"
    end

    # Categorize primes
    for n in SELECTED_PRIMES
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
                CUMULATIVE_PRODUCTS, SELECTED_PRIMES)
end
