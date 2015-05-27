from parametrized_set import ParametrizedSet
from sympy import integrate, Expr, Number, NumberSymbol, S

#FIXME 0-measure
#FIXME Dirac measure

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

    def __add__(self, other):
        '''Product of two measures is a new ProductMeasure.'''
        assert isinstance(other, (Measure, ProductMeasure))
        
        if isinstance(other, ProductMeasure):
            return ProductMeasure([self] + other.measures)
        else:
            return ProductMeasure([self, other])


class ProductMeasure(Measure):
    '''Product of measures.'''

    def __init__(self, measures):
        '''Initialize from the list of measures.'''
        assert isinstance(measures, list)

        # Check domain compatibility. this does not mean that integration won't 
        # blow up
        gdim = measures[0].domain.gdim
        tdim = measures[0].domain.tdim
        assert all(measure.domain.tdim == tdim and measure.domain.gdim == gdim
                   for measure in measures[1:]),\
                           'Cannot sum measures of different tdim and gdim'
        self.measures = measures

        # No explicit domain
        Measure.__init__(self, None)

    def __call__(self, integrand):
        '''Call makes no sense with None domain.'''
        raise NotImplementedError('No __call__ for product measure')

    def __rmul__(self, integrand):
        '''Integrate with individual measures.'''
        measures = self.measures
        ans = integrand*measures[0]
        # FIXME if we add += to vectors/tensors this can be simplified
        for measure in measures[1:]:
            ans = ans + integrand*measure
        return ans

    def __add__(self, other):
        '''Combine measures.'''
        assert isinstance(other, (Measure, ProductMeasure))

        if isinstance(other, ProductMeasure):
            return ProductMeasure(self.measures + other.measures)
        else:
            return ProductMeasure(self.measures + [other])

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
