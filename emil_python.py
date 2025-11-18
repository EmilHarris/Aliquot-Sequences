# Imports
import math as m
import numpy as np

# Set default values for K, I, and the maximum n value
K_norm = 30

I_norm = 1e9

n_max = 20000


# Prime sieve is an efficient way to find all the primes up to (and including) n
def prime_sieve(maximum=n_max):

    # There are no primes less than 2
    if maximum < 2:
        return None

    # Will create a list of True or False values corresponding to whether the index is a prime number
    # e.g. 3 is a prime number so the value in the 3rd index (4th position) will be true by the end of the function
    prime_mask = [True for _ in range(maximum+1)]

    # 0 and 1 are not prime
    prime_mask[0], prime_mask[1] = False, False

    # Only need to check up to root n
    for p in range(2, int(m.sqrt(maximum)) + 1):
        if prime_mask[p]:
            # go over multiples of each prime number, which will not be prime, and so set value in prime_mask to False
            for i in range(p*p, maximum+1, p):
                prime_mask[i] = False

    # change prime_mask into a numpy array so we can use a boolean mask
    prime_mask = np.array(prime_mask)

    # generate the array of numbers up to n and apply boolean mask to create array of prime numbers
    return np.array([i for i in range(maximum+1)])[prime_mask]


# PRIMES calculated using prime sieve up to the largest n we are going up to, divided by 2
# We can stop here as n cannot have a prime factor larger than n/2
# this can then be used to calculate s(n) for ALL the ns up to n_max, instead of prime_sieve being called for every n
def get_primes(n=n_max):

    try:
        with open('primes.txt', 'r') as f:

            lines = f.readlines()
            n_cached = lines[0]
            primes = [int(p) for p in lines[1:]]

        if n_cached < n:
            primes.append(prime_sieve())

    except FileNotFoundError:
        primes = prime_sieve(int(n/2) + 1)



    return primes


PRIMES = get_primes()

# Using the PRIMES list get the prime factors of a umber n
def get_prime_factors(n):

    primes = PRIMES

    # create an empty dictionary of prime factors
    # we will store this is prime : exponent pairs
    # e.g. n=12 will have prime_factors = {2 : 2, 3 : 1} as 12 = 2^2 x 3
    prime_factors = {}

    # go through each prime to check if it's a divisor of n
    for prime in primes:

        # Only need to check primes up to root n
        # There may be 1 prime bigger than root n,we will deal with this at the end
        if prime**2 > n:
            break

        # check if prime divides n and add to dictionary if it does
        if n % prime == 0:
            prime_factors[int(prime)] = 1

            # now we set n to n / prime so we dont count the same prime factor more than once
            n = n / prime

            # check if this new n is still divisible by the same prime.
            # For as long as it works keep increasing the exponent value assigned to the prime, and dividing n
            while n % prime == 0:
                prime_factors[int(prime)] += 1
                n = n / prime

        # if n reaches 1 we can stop early as we know we have found the prime factorisation of n
        if n == 1:
            break

    # now we check for the potential prime factor bigger than root n.
    # If we have checked all the primes up to root n and still haven't found all the prime factors,
    # we must be left with the last prime factor as n
    if n > 1:
        prime_factors[int(n)] = 1

    # return our completed list of prime factors
    return prime_factors


# Calculate s(n)
def s(n):

    # 1 has no proper divisors
    if n == 1:
        return 0

    # the following uses the formula: (insert in latex later)
    # to calculate the sum of all the factors of n using only its prime factors
    total = 1

    primes = get_prime_factors(n).items()

    for prime, power in primes:
        total *= (prime**(power+1) - 1)/(prime - 1)

    # return s(n)
    return total - n


# Calculate the aliquot sequence of n
def aliq_seq(n, K=K_norm, I=I_norm):

    # n will be the first element of its aliquot sequence
    seq = [n]

    # make an empty set to record all the s(n)s in order to detect loops
    seen = set()

    curr_n = n
    # loop max k times to compute and check each s(n)
    # return always ends the function so no more code will run after one of the checks is 'failed'
    for k in range(K - 1):

        # check if sequence is terminated
        if curr_n == 0:
            return seq, 'terminated'

        # check if number is bigger than specified i
        if curr_n >= I:
            return seq, 'I reached'

        # check if sequence has looped
        if curr_n in seen:
            return seq, 'looped'

        # add n to sequence of seen numbers
        # has to be done after checks otherwise it will always think it loops
        seen.add(curr_n)

        # calculate next term of sequence if n has 'passed' all tests
        curr_n = s(curr_n)
        seq.append(int(curr_n))

    # If loop finished and has computed all k s(n), the sequence is finished (up to where we have decided to stop)
    # so returns sequence and the status that k has been reached
    return seq, 'K reached'


# make dictionary to count the number of types of each sequence
counts = {'terminated': 0, 'looped': 0, 'I reached': 0, 'K reached': 0}

# generate sequences up to 20000 (or whatever n_max is) and update counts dictionary
for i in range(1, n_max + 1):
    status = aliq_seq(i)[1]
    counts[status] += 1

# output counts (for testing)
print(counts)
