# Imports
import math as m
import numpy as np
import time as t
import json

# Set default values for K, I, and the maximum n value
K_norm = 60

I_norm = 1e11

n_max = 100000

json_max = 30000

start_time = t.perf_counter()


def extend_prime_mask(old_mask, old_max, new_max):
    prime_mask = old_mask + [True for _ in range(old_max + 1, new_max + 1)]

    prime_mask[0], prime_mask[1] = False, False

    for p in range(2, int(m.sqrt(new_max)) + 1):
        if prime_mask[p]:

            for i in range(p**2, new_max+1, p):
                prime_mask[i] = False

    return prime_mask


# PRIMES calculated using prime sieve up to the largest n we are going up to, divided by 2
# We can stop here as n cannot have a prime factor larger than n/2
# this can then be used to calculate s(n) for ALL the ns up to n_max, instead of prime_sieve being called for every n
def get_prime_mask(n=n_max):

    try:
        f = open('prime_mask.txt', 'r')

        lines = f.readlines()
        prime_mask = [bool(x) for x in lines]
        maximum = len(prime_mask) - 1

        f.close()

    except FileNotFoundError:
        prime_mask = []
        maximum = -1

    if maximum < n:
        prime_mask = extend_prime_mask(prime_mask, maximum, n)

    prime_str = [str(x) for x in prime_mask]

    with open('primes.txt', 'w') as f:
        f.write('\n'.join(prime_str))

    return prime_mask[:n+1]


nums = np.array([i for i in range(n_max + 1)])
PRIMES = nums[get_prime_mask()]
print(f'primes found after {t.perf_counter() - start_time}')


def convert_keys(d):
    new = {}
    for k, v in d.items():
        k = float(k)
        new[k] = convert_keys(v) if isinstance(v, dict) else v
    return new


try:
    with open('factorisations.json', 'r') as f:
        factorisations = convert_keys(json.load(f))

except FileNotFoundError:
    factorisations = {}


# Using the PRIMES list get the prime factors of a umber n
def get_prime_factors(n):

    global factorisations

    if n in factorisations:
        return factorisations[n]

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

            if n in factorisations:
                facts = factorisations[n].copy()

                for key, power in prime_factors.items():
                    try:
                        facts[key] += power

                    except KeyError:
                        facts[key] = power

                return facts

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

    global factorisations

    # 1 has no proper divisors
    if n == 1:
        return 0

    # the following uses the formula: (insert in latex later)
    # to calculate the sum of all the factors of n using only its prime factors
    total = 1

    primes_dict = get_prime_factors(n)
    primes = primes_dict.items()
    factorisations[n] = primes_dict

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

print(f'aliquot sequences found after {t.perf_counter() - start_time}')
print(counts)

with open('factorisations.json', 'w') as f:
    json_factorisations = {}

    for key in factorisations.keys():
        if key <= json_max:
            json_factorisations[int(key)] = factorisations[key]
    factorisations_obj = json.dumps(factorisations)
    f.writelines(factorisations_obj)

# output counts (for testing)
print(f'finished after {t.perf_counter() - start_time}')

