from scipy.stats import chi2

def confidence_interval_probability_bernoulli_model_independent_examinations(n, m, beta):
    def R0(n, beta):
        return n*(1-pow(1-beta, 1/n))
    def R1(n, m, beta):
        x = chi2.ppf(1-beta, 2*m)
        # print(f"x(1-b) = {x}, df = {2*m}; q = {1-beta}")
        res = (m*(2*n-m+1+1/2*x))/(n*x)
        return res
    def R2(n, m, beta):
        x = chi2.ppf(beta, 2*(m+1))
        # print(f"x(b) = {x}, df = {2*(m+1)}; q =  {beta}")
        return (m*(2*n-m+1/2*x))/(n*x)
    if m == 0:
        p1 = 0
        p2 = R0(n, beta)/n
    else:
        p1 = m/(n*R1(n, m, beta))
        p2 = m/(n*R2(n, m, beta))
    return (p1, p2)

