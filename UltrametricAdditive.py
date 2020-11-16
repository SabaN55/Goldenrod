from compsci260lib import *


def solve_ultrametric_additive():
    # Distance metrics for table 1 and table 2
    dist_1 = {"1,2": 0.3, "1,3": 0.7, "1,4": 0.9,
              "2,3": 0.6, "2,4": 0.8,
              "3,4": 0.6}

    dist_2 = {"1,2": 0.8, "1,3": 0.4, "1,4": 0.6, "1,5": 0.8,
              "2,3": 0.8, "2,4": 0.8, "2,5": 0.4,
              "3,4": 0.6, "3,5": 0.8,
              "4,5": 0.8}

    # Check if dist_1 and dist_2 are ultrametric and additive by
    # calling is_ultrametric and is_additive with the default
    # threshold value (1e-4).
    #
    print(is_ultrametric(dist_1, threshold=1e-4))
    print(is_ultrametric(dist_2, threshold=1e-4))
    print(is_additive(dist_1, threshold=1e-4))
    print(is_additive(dist_2, threshold=1e-4))
    #    

    # Construct the ATPA synthase distance metric table
    atpa_table = {"1,2": 0.5, "1,3": 0.5, "1,4": 0.1, "1,5": 0.4, "1,6": 0.4,
                  "2,3": 0.3, "2,4": 0.5, "2,5": 0.5, "2,6": 0.5,
                  "3,4": 0.5, "3,5": 0.5, "3,6": 0.5,
                  "4,5": 0.4, "4,6": 0.4,
                  "5,6": 0.3}

    # Determine if the ATPA synthase distance metrics
    # are ultrametric and additive using the default
    # threshold value (1e-4).
    #
    print(is_ultrametric(atpa_table, threshold=1e-4))
    print(is_additive(atpa_table, threshold=1e-4))
    #    


def is_ultrametric(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are ultrametric.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal
    Returns:
        (bool) True if the given distance metric is an ultrametric,
    False otherwise."""

    #
    # First we start by finding the total number of points
    num_points = num_points_finder(dist)
    p_range = range(1, num_points + 1)

    tocomp = []
    # Now we iterate over all these points and find every unique 3-way connection between them
    for i in p_range:
        for j in p_range:
            if i != j:
                for k in p_range:
                    if k != i and k != j:
                        toadd = sorted([i, j, k])
                        if toadd not in tocomp:
                            tocomp.append(sorted([i, j, k]))
    # print(tocomp)

    # Here we iterare through the list of connections and check to make sure that each
    # of them follows the rules of an ultrametric tree, if they do not, a False is returned.
    for [i, j, k] in tocomp:
        triplet = [i, j, k]
        case = 0
        key1 = f"{i},{j}"
        key2 = f"{i},{k}"
        key3 = f"{j},{k}"

        d1 = dist.get(key1)
        d2 = dist.get(key2)
        d3 = dist.get(key3)

        ascending = sorted([d1, d2, d3])

        if is_almost_equal(ascending[1], ascending[2], threshold):
            if ascending[0] >= ascending[1]:
                print(f"The triplet {triplet} fails the criteria.")
                return False
        else:
            print(f"The triplet {triplet} fails the criteria.")
            return False
    return True


def is_additive(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are additive.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal

    Returns:
        (bool) Return True if the given distance metric is additive, 
        False otherwise."""

    #
    num_points = num_points_finder(dist)
    p_range = range(1, num_points + 1)

    tocomp = []
    # Now we iterate over all these points and find every unique 3-way connection between them
    for i in p_range:
        for j in p_range:
            if i != j:
                for k in p_range:
                    if k != i and k != j:
                        for l in p_range:
                            if l != i and l != j and l != k:
                                toadd = sorted([i, j, k, l])
                                if toadd not in tocomp:
                                    tocomp.append(sorted([i, j, k, l]))
    # print(tocomp)

    # Here we iterare through the list of connections and check to make sure that each
    # of them follows the rules of an ultrametric tree, if they do not, a False is returned.
    for [i, j, k, l] in tocomp:
        quartet = [i, j, k, l]
        case = 0
        key1 = f"{i},{j}"
        key2 = f"{i},{k}"
        key3 = f"{i},{l}"
        key4 = f"{j},{k}"
        key5 = f"{j},{l}"
        key6 = f"{k},{l}"

        d1 = dist.get(key1)  # ij
        d2 = dist.get(key2)  # ik
        d3 = dist.get(key3)  # il
        d4 = dist.get(key4)  # jk
        d5 = dist.get(key5)  # jl
        d6 = dist.get(key6)  # kl

        comp1 = d1 + d6  # ij+kl
        comp2 = d2 + d5  # ik+jl
        comp3 = d3 + d4  # il+jk

        ascending = sorted([comp1, comp2, comp3])

        if is_almost_equal(ascending[1], ascending[2], threshold):
            if ascending[0] >= ascending[1]:
                print(f"The quartet {quartet} fails the criteria.")
                return False
        else:
            print(f"The quartet {quartet} fails the criteria.")
            return False
    return True


def is_almost_equal(num_1, num_2, threshold):
    """
    Return true if the difference between the two parameters is negligible
    enough that the parameters can be considered equal.
    """
    return abs(num_1 - num_2) <= threshold


def num_points_finder(dict):
    '''
    This helper function allows me to calculate the number of points that I will have to
    compare.
    '''
    pairs = dict.keys()
    max = 0
    for item in pairs:
        tocomp = item.split(",")
        bigboy = int(tocomp[1])
        if bigboy > max:
            max = bigboy
    return max


if __name__ == '__main__':
    solve_ultrametric_additive()
