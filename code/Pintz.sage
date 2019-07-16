print("Siegel zeroes for primitive real characters.")

even = 'even'
odd = 'odd'

var('t')

def character_sum(q0, q1, parity):      # Lapkova 2018
    """Calculate upper bound for sum of Dirichlet characters of modulus q0 <= q <= q1."""

    if parity == 'even':
        number = (2/pi**2 + 0.9467/log(q0) + 1.668/sqrt(q0)/log(q0))
    if parity == 'odd':
        number = (0.5/pi + 0.8294/log(q0) + 1.0285/sqrt(q0)/log(q0))
    return number*sqrt(q1)*log(q1)

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

    alpha = 1 - tau

    def h6(d):
        return log(d)/d**(1-tau)
    def h7(d):
        return 1/d**(1-tau)
    def h8(d):
        return 1/d

    # d/dx H[i](x) = h[i](x)

    def H6(d):
        return d**tau/tau*(log(d) - 1/tau)
    def H7(d):
        return d**tau/tau
    def H8(d):
        return log(d)

    K1 = abs(C1)
    K2 = abs(C2)
    K3 = abs(x^tau/tau*(1/tau - log(x)))

    # Spitting sum
    D1 = x**tau*log(x)/2/x     # lead term of W
    D2 = A/tau**2              # lead term of upper_sum
    z = sqrt(D2/D1)            # minimize D1*z + D2/z
    z = min(z, x/2)
    z = floor(z)

    L_error1, L_error2, L_error3 = +Infinity, +Infinity, +Infinity

    if parity == 'even':        # Louboutin 2001
        N0 = floor(sqrt(q0)) - 1
        N1 = floor(sqrt(q1)) - 1
        L_error1 = (H6(N1 + z) - H6(z) + N1*(N1 + 1)/2*h6(N0 + z + 2)
                    - N0*(N0+3)/2*h6(N1 + z + 1)
                    + sqrt(q1)/2*((N1 + 1)*h6(N0 + z + 1) - N0*h6(N1 + z + 2)))
        L_error2 = (H7(N1 + z) - H7(z) + N1*(N1 + 1)/2*h7(N0 + z + 2)
                    - N0*(N0+3)/2*h7(N1 + z + 1)
                    + sqrt(q1)/2*((N1 + 1)*h7(N0 + z + 1) - N0*h7(N1 + z + 2)))
        L_error3 = (H8(N1 + z) - H8(z) + N1*(N1 + 1)/2*h8(N0 + z + 2)
                    - N0*(N0+3)/2*h8(N1 + z + 1)
                    + sqrt(q1)/2*((N1 + 1)*h8(N0 + z + 1) - N0*h8(N1 + z + 2)))

    L_error1 = min(L_error1, H6(A + z) - H6(z-1))
    L_error2 = min(L_error1, H7(A + z) - H7(z-1))
    L_error3 = min(L_error1, H8(A + z) - H8(z-1))

    W = (x**tau*log(x)/2*z/x
         + x**tau*(1 + alpha*log(x))/24*z/x*(z/x + 1/x)
         + x**tau*(alpha*(alpha + 1)*log(x) + (3*alpha**2 + 6*alpha + 2)/(alpha + 2))
         /216/sqrt(3)*z/x*(z/x + 1/x)*(2*z/x + 1/x))

    # W1 = x**tau*log(x)/2*z/x
    # W2 = x**tau*(1 + alpha*log(x))/24*z/x*(z/x + 1/x)
    # W3 = (x**tau*(alpha*(alpha + 1)*log(x) + (3*alpha**2 + 6*alpha + 2)/(alpha + 2))
    #       /216/sqrt(3)*z/x*(z/x + 1/x)*(2*z/x + 1/x))

    upper_sum = A/z**(1-tau)*(log(z) + (x/z)**tau/tau*(log(x) - 1/tau)
                              - (log(z) - 1/tau)/tau)

    # print(float(W1), float(W2), float(W3), float(upper_sum))

    # print(three_six.n(), three_seven.n(), three_eight.n(), W.n(), upper_sum.n())

    if 4e5 <= q0 and q1 <= 1e7:
        Bennet = 79.2
    else:
        Bennet = 12.52    # Bordignon

    L1 = (1/tau - log(z))*Bennet*z^tau/tau/sqrt(q1)
    lower_bound = max(-2*zetaderiv(1, 2-2*tau)
                      -2*(1 + (0.5 - tau)*log(z))/z^(0.5 - tau)/(1 - 2*tau)^2,
                      log(4)/4)

    number = (L_error1 + L_error2 + L_error3 + W + upper_sum
              - L1 - lower_bound)
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
    x_max = min(x1, log(q0^(1/c), 10))    # upper bound from (10)
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
