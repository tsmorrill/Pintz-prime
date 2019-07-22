print("Siegel zeroes for primitive real characters of prime modulus.")

even = 'even'
odd = 'odd'

var('t')

def character_sum(q0, q1, parity):      # Lapkova 2018
    """Calculate upper bound for sum of Dirichlet characters of modulus q0 <= q
    <= q1. This should be multiplied by the square root of the length of the sum
    .
    """

    even_table = [1.83, 1.74, 1.68, 1.63, 1.59, 1.55, 1.52, 1.49, 1.47, 1.45,
                  1.43, 1.41, 1.39, 1.38, 1.12, 1.02]
    odd_table = [1.79, 1.70, 1.64, 1.60, 1.56, 1.53, 1.50, 1.47, 1.45, 1.43,
                 1.41, 1.39, 1.38, 1.36, 1.10, 1.00]

    index = int(log(q0, 10) - 7)
    index = min(index, len(even_table)-1)

    if parity == 'even':
        constant = even_table[index]
    if parity == 'odd':
        constant = odd_table[index]

    return constant*q1**(3/16)*sqrt(log(q1))

def constants(alpha):
    """Calculate constants related to Euler-Maclaurin summation for n**-alpha
    and n**-alpha*log(n)."""

    def B3(t):                # Third Bernoulli polynomial
        t = frac(t)
        number = t**3 - 1.5*t**2 + 0.5*t
        return number

    integral1 = -(alpha*(alpha + 1)*(alpha + 2)/6
                  *numerical_integral(B3(t)/t**(alpha + 3), 1, Infinity)[0])
    integral2 = 1/6*numerical_integral(B3(t)/t**(alpha + 3)
                   *(alpha*(alpha+1)*(alpha+2)*log(t)
                     + 3*alpha**2 + 6*alpha + 2), 1, Infinity)[0]

    C1 = (0.5 + alpha/12 - 1/(1-alpha) + integral1)
    C2 = (-1/12 + 1/(1 - alpha)**2 + integral2)

    # print(C1, C2)
    return(C1, C2)

def error(q0, q1, A, tau, x, C1, C2, parity):
    """Calculate F on the interval [q0, q1] for precomputed A, tau, C1, C2."""

    upper_bound = 6*x**tau*log(x)*sqrt(A)/t**2/x

    if 4e5 <= q0 and q1 <= 1e7:
        Bennet = 79.2
    else:
        Bennet = 12.52    # Bordignon

    L1 = (1/tau - log(z))*Bennet*z^tau/tau/sqrt(q1)
    lower_bound = max(-2*zetaderiv(1, 2-2*tau)
                      -2*(1 + (0.5 - tau)*log(z))/z^(0.5 - tau)/(1 - 2*tau)^2,
                      log(4)/4)
    lower_bound_classic = log(4)/4

    number = (upper_bound -  lower_bound_classic)
    return number

def F(c, q0, q1, x, parity):
    """Calculate an upper bound for F on the interval [q0, q1] for fixed c and x."""

    A = character_sum(q0, q1, parity)
    tau = c/log(q1)
    alpha = 1 - tau
    C1, C2 = constants(alpha)
    number = error(q0, q1, A, tau, x, C1, C2, parity)

    return (float(number))

def search(c, q0, q1, x0, x1, step, parity):
    """Find an x on the interval [x0, x1] which minimizes F for c and q in [q0, q01],
    searching in increments of 10**step. Then, return the tuple (bool(F < 0), log(x, 10)).
    Returns None if the specified range and step size have no test values.
    """

    tau = c/log(q1)
    A = character_sum(q0, q1, parity)

    x_min = max(x0, 2, log((exp(1/(2*tau - 1)) + 1)**2, 10))    # lower bound from Lemma 4; also need log(x) >0
    x_max = min(x1, log(q0, 10)/4)    # upper bound from (10)
    x_range = [x_min + i*step for i in range(floor((x_max-x_min)/step) + 1)]
    values = []

    if not x_range:
        return (False, 10)   # x value should be ignored here

    alpha = 1 - tau
    C1, C2 = constants(alpha)

    for log_x in x_range:
        x = 10^log_x
        number = error(q0, q1, A, tau, x, C1, C2, parity)
        values.append((float(number), log_x))
    F, x = (float(i) for i in min(values))
    output = (F < 0)

    return(output, x)

def best_c(q0, q1, significant_figures=2):
    """Calculate the maximal c so that F < 0 on the interval [q0, q1]."""

    parity = 'even'
    print('Sigel zeroes for even characters.')
    print('')

    c_even, c_step = 0, 1.0
    record_significant_figures, current_figures = False, 0
    done = False

    while not done:
        print('Trying increments of {}.'.format(c_step))
        c, it_works = c_even + c_step, True
        while it_works:
            x0, x1, step = 0, log(q0^(1/c), 10), 1
            for i in range(3):
                result, x = search(c, q0, q1, x0, x1, step, parity)
                x0, x1, step = max(0, x - step), min(x1, x + step), step/10
                if result:
                    it_works = True
                    record_significant_figures = True
                    x_even, c_even = x, c
                    print("q in [{}, {}], c >= {}.".format(q0, q1, c_even))
                    c += c_step
                    if c > 4:
                        done = True
                        print('ERROR: c > 4.')
                    break
                else:
                    it_works = False
        if record_significant_figures == True:
            current_figures += 1
            if current_figures == significant_figures:
                done = True
        c_step /= 10
        print('')

    parity = 'odd'
    print('Sigel zeroes for odd characters.')
    print('')

    c_odd, c_step = 0, 1.0
    record_significant_figures, current_figures = False, 0
    done = False

    while not done:
        print('Trying increments of {}.'.format(c_step))
        c, it_works = c_odd + c_step, True
        while it_works:
            x0, x1, step = 0, log(q0^(1/c), 10), 1
            for i in range(3):
                result, x = search(c, q0, q1, x0, x1, step, parity)
                x0, x1, step = max(0, x - step), min(x1, x + step), step/10
                if result:
                    it_works = True
                    record_significant_figures = True
                    x_odd, c_odd = x, c
                    print("q in [{}, {}], c >= {}.".format(q0, q1, c_odd))
                    c += c_step
                    break
                else:
                    it_works = False
        if record_significant_figures == True:
            current_figures += 1
            if current_figures == significant_figures:
                done = True
        c_step /= 10
        print('')

    q0_magnitude = int(log(q0, 10))
    q0_lead = round(q0/10**q0_magnitude, 2)
    q1_magnitude = int(log(q1, 10))
    q1_lead = round(q1/10**q1_magnitude, 2)

    c_even, c_odd = round(c_even, 5), round(c_odd, 5)

    TeX_string = "${{{}}} \cdot 10^{{{}}}$ & ${{{}}} \cdot 10^{{{}}}$ & \\num{{{}}} & $10^{{{}}}$  & \\num{{{}}} & $10^{{{}}}$  \\\\".format(q0_lead, q0_magnitude, q1_lead, q1_magnitude, c_even, x_even, c_odd, x_odd)
    even_string = 'F({}, {}, {}, 10^{}, {})'.format(c_even, q0, q1, x_even, even)
    odd_string = 'F({}, {}, {}, 10^{}, {})'.format(c_odd, q0, q1, x_odd, odd)

    print('')
    print('For even chi, c = {}.'.format(c_even))
    print('For odd chi,  c = {}.'.format(c_odd))
    print('')
    return(TeX_string, even_string, odd_string)


def cq_table(q_list, significant_figures=4):
    """Generate a table of c values corresponding to q_list formatted for LaTeX, then generate Sage
    commands to verify these calculations.
    """

    TeX_list, even_list, odd_list = [], [], []
    for q0, q1 in zip(q_list[:-1], q_list[1:]):
        (TeX_string, even_string, odd_string) = best_c(q0, q1, significant_figures=significant_figures)
        print('====')
        TeX_list.append(TeX_string)
        even_list.append(even_string)
        odd_list.append(odd_string)

    print('')
    print('LaTeX Table:')
    for TeX_string in TeX_list:
        print(TeX_string)
    print('')
    print('Verification for even characters:')
    for even_string in even_list:
            print(even_string)
    print('')
    print('Verification for odd characters:')
    for odd_string in odd_list:
        print(odd_string)

def subdivide(list):
    new_list = []
    for q0, q1 in zip(list, list[1:]):
        number = sqrt(q0*q1)
        number_magnitude = int(log(number, 10))
        number_lead = round(number/10**number_magnitude, 1)
        new_list.append(number_lead*10**number_magnitude)
    new_list += list
    new_list.sort()
    return new_list
