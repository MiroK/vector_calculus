from parametrized_set import ParametrizedSet
from sympy import integrate, Expr, Number, NumberSymbol, S


class Measure(object):
    '''Integral over domain describing points in Cartesian coordinate system.'''

    def __init__(self, domain):
        self.domain = domain

    def __call__(self, integrand):
        '''
        Integrate scalar integrand with the measure. This is a working horse for 
        specialized classes. There is no Jacobian!
        '''
        if isinstance(integrand, (int, float)):
            integrand = S(integrand)
        
        assert isinstance(integrand, (Expr, Number, NumberSymbol))

        # Substitute
        f = self.domain.substitute(integrand)

        # Integrate over parameter domain
        ans = f
        for var, bounds in self.domain.items():
            ans = integrate(ans, (var, bounds[0], bounds[1]))

        return ans
        

# -----------------------------------------------------------------------------
        

if __name__ == '__main__':
    from parametrized_set import Interval
    from sympy import Symbol

    
    domain = Interval(-1, 1)
    dl = Measure(domain)
    
    x = Symbol('x')
    f = 1 + x**2

    print dl(f)
    print float(integrate(f, (x, -1, 1)))
