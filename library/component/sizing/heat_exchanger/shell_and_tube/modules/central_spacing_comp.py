#%% CENTRAL SPAC RELATED FUN

import numpy as np
import random

# Function to find divisors between bounds for floating-point numbers
def find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound, tolerance=1e-9):
    """ Finds floating-point divisors of num_to_divide within a specific range. """
    valid_divisors = []
    
    # Loop through potential divisors in small increments between lower and upper bounds
    potential_divisors = np.arange(round(lower_bound, 2), round(upper_bound, 2), 0.01)  # Adjust step size as needed
    
    print(round(lower_bound, 2))
    print(round(upper_bound, 2))

    for divisor in potential_divisors:
        divisor_rounded = round(divisor,2)
        if np.mod(num_to_divide,divisor_rounded) == 0:
            valid_divisors.append(divisor_rounded)
    
    print(f"valid_divisors : {valid_divisors}")

    return valid_divisors

# Function to generate a random divisor between two bounds
def random_divisor_between_bounds(num_to_divide, lower_bound, upper_bound):
    valid_divisors = find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound)
    
    if valid_divisors:
        # Choose a random divisor from the valid divisors list
        return random.choice(valid_divisors)
    else:
        return None  # No divisors found in the given range
