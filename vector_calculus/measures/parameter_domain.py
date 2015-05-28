from collections import OrderedDict
from sympy import Symbol, Number, NumberSymbol


class ParameterDomain(object):
    '''Specify domain of parametrized sets.'''

    def __init__(self, *args):
        '''Domain is build from tuple of pairs.'''
        self._domain = OrderedDict({})
        # Vars in terms of which the domain is defines
        free_vars = set([])
        for arg in args:
            var = arg[0]
            # Parameters must be unique symbols
            assert isinstance(var, Symbol), 'Parameter must by Symbol'
            assert var not in free_vars,\
                    'Variable %s already in free variables %s' % (str(var), free_vars)
           
            # Limits can be only functions of old parameters
            assert len(arg[1]) == 2, 'Bounds must have length 2'
            # Turn the arguments to sympy
            limits = map(lambda v: Number(v) if isinstance(v, (float, int)) else v,
                         arg[1])
            # Check that the limits are not defined in terms of current
            # parameter
            assert all(var not in l.atoms() for l in limits), 'No implicit definitions'
            # Check that the limts are either independent of functions of free
            # parameters
            for l in limits:
                # If the conversion fails we have an Expr to check
                try:
                    float(l)
                except TypeError:
                    assert len(l.atoms().intersection(free_vars)) > 0, \
                        'Bounds must be defined in terms of existing parameters'

            # Add the new parameter if everything okay
            free_vars.add(var)

            # Also to dictionary
            self._domain[var] = limits

            # Unknowns
            self._parameters = set(self._domain.keys())

    def __getitem__(self, var):
        '''Bounds for the var parameter.'''
        return self._domain[var]

    def __iter__(self):
        '''Iterate in reversed order (order of integration).'''
        return self.items()

    def __len__(self):
        '''Something like a topological dimension.'''
        return len(self._domain)

    def items(self):
        '''Reversed items iterator.'''
        return reversed(self._domain.items())

    @property
    def parameters(self):
        '''Parameters that define the domain.'''
        return self._parameters

# -----------------------------------------------------------------------------


if __name__ == '__main__':
    from sympy import symbols

    s, t, r = symbols('s, t, r')

    # Simplex like
    line_pdomain = ParameterDomain((s, (0, 1)))
    triangle_pdomain = ParameterDomain((s, (0, 1)), (t, (0, 1-s)))
    tetrahedron_pdomain = ParameterDomain((s, (0, 1)),
                                          (t, (0, 1-s)),
                                          (r, (0, 1-s-t)))

    # Cartesian product like
    interval_pdomain = ParameterDomain((s, (-1, 1)))
    rectangle_pdomain = ParameterDomain((s, (-1, 1)), (t, (-1, 1)))
    box_pdomain = ParameterDomain((s, (-1, 1)), (t, (-1, 1)), (r, (-1, 1)))
