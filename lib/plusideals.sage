# -*- coding: utf-8 -*-

def ab_slope(a, b, cloud):
    return QQ(cloud[b][1]-cloud[a][1])/QQ(cloud[b][0]-cloud[a][0])

def lower_convex_hull(cloud):
    """
    Computes the lower convex hull of a cloud of points.
    """
    slopes = [ [0, 1, ab_slope(0, 1, cloud)] ]
    b = 2
    while b < len(cloud):
        for i in [0..len(slopes)-1]:
            a = slopes[i][0]
            slope = ab_slope(a, b, cloud)
            if slope <= slopes[i][2]:
                slopes = slopes[0:i]
                slopes.append([a, b, slope])
                break
        if slopes[-1][1] != b:
            a = b - 1
            slope = ab_slope(a, b, cloud)
            slopes.append([a, b, slope])
        b += 1

    slopes = [ [slope, cloud[a], cloud[b]] for a, b, slope in slopes ]
    return slopes


class MontesType(object):

    def __init__(self, K, p, varphi, omega, initial=True):
        self.parent = K
        self.pol = K.defining_polynomial()
        self.prime = p
        self.varphi = varphi
        self.levels = [ ]
        if initial is True:
            self.add_initial_level(varphi, omega, p)
    
    def copy(self):
        copy = MontesType(self.parent, self.prime, self.varphi, self.levels[0].omega, initial=False)

        for level in self.levels:
            copy.levels.append(level.copy())
        
        return copy

    def add_initial_level(self, varphi, omega, p):
        phi = self.pol.parent()(varphi)

        new_level = MontesTypeLevel(phi, omega, p)

        new_level.Fq = FiniteField(p^varphi.degree(), 'z0', varphi)
        new_level.z = new_level.Fq.0
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y')

        self.levels.append(new_level)

    def add_new_level(self, phi, omega, u):
        s = len(self.levels)
        lvl_s = self.levels[-1]
        new_level = MontesTypeLevel(phi, omega, self.prime)

        new_level.V = lvl_s.e^u
        new_level.prod_e = lvl_s.prod_e * lvl_s.e
        new_level.prod_f = lvl_s.prod_f * lvl_s.f
        
        # FIXME: p^f, f = what??
        new_level.Fq = FiniteField(self.prime^new_level.prod_f, 'z'+str(s), lvl_s.res_pol)
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y'+str(s))

        if lvl_s.f > 1:
            new_level.z = new_level.Fq.0
        else:
            new_level.z = lvl_s.res_pol.coefficients()[0]

        print "Fq: %s, Fpy: %s, z: %s" % (str(new_level.Fq), str(new_level.Fqy), str(new_level.z),)
        
        self.levels.append(new_level)
        return new_level

    def __unicode__(self):
        return '(%s)' % (', '.join([unicode(l) for l in self.levels]))

    def __repr__(self):
        return self.__unicode__()

    def rth_level(self):
        return self.levels[-1]


    #### Begin serious Montes methods ####

    ## Phi-Newton Polygons ##
    def phi_newton_polygon(self, i, phiadic):
        assert i <= len(self.levels)

        n = 0
        sides = []
        cloud = []
        side_devs = []
        all_devs = []
        for k in [0..len(phiadic)-1]:
            val = 0
            dev = [ ]

            val, dev = self.value(i, phiadic[k])
            #print "N_%d : v_%d(%s phi_%d) = %s + %d" % (i, i, phiadic[k], i, str(val), n)
            if abs(val) != Infinity:
                cloud.append([k, val + n])
                all_devs.append(dev)
            n += self.levels[i-1].V

        sides = lower_convex_hull(cloud)
        abscissas = [ p[0] for p in cloud ]    

        for side in sides:
            height = ZZ(side[1][1])
            dev = [ ]

            for j in range(side[1][0], side[2][0] + side[0].denominator(), side[0].denominator()):
                try:
                    position = abscissas.index(j)
                except ValueError:
                   position = -1 
                
                if position > -1 and cloud[position][1] == height:
                    dev.append(all_devs[position])
                elif i == 1:
                    dev.append(0)
                else:
                    dev.append([])
                height += side[0].numerator()
            dev.append(list(side[1]))
            side_devs.append(dev)    

        if len(sides) == 0:
            raise NotImplemented, "Case for a zero side Newton Polygon has not been implemented"


        return sides, side_devs

    ## Value of polynomial a(x) in ZZ[x] ##
    def value(self, i, a):
        assert i <= len(self.levels)

        val = +Infinity
        if i == 1:
            if a != 0:
                val = min([valuation(c, self.prime) for c in a.coefficients()])
            devs = a
        else:
            if a == 0:
                devs = [ ]
            else:
                raise NotImplemented, "Valuation for i > 1 not implemented."

        return val, devs

    ## Residual Polynomial ##
    def residual_polynomial(self, i, side_devs):
        """
        Creates and returns the i-th residual polynomial of a polynomial (f).
        side_devs is a list of the multiadic expansions of the coefficients of
        f whose attached points in N_r(f) lie on teh side S of slope
        -levels[i].slope. The last element of side_devs is [s, u], where s,u
        are the coordinates of the left end-point on S.
        """

        assert i <= len(self.levels)

        lvl_i = self.levels[i-1]
        height = side_devs[-1][1]
        res_coeffs = [ ]

        for dev in side_devs[:-1]:
            if (i == 1 and dev == 0) or (i > 1 and len(dev) == 0): 
                res_coeffs.append(lvl_i.Fq(0))
            else:
                if i == 1:
                    c = ZZ[x](ZZ(dev) // self.prime**height)
                    res_coeffs.append(c(lvl_i.z))
                    j = side_devs.index(dev)
                    print "%d. order Res.Pol. the easy way c_%d = (%d / %d) = %d (z_%d = %d)" % (i, j, dev, self.prime**height, dev // self.prime**height, i-1, lvl_i.z)
                else:
                    raise NotImplemented, "Level i > 1 res pols are not yet implemented."
            height = height - lvl_i.h

        return lvl_i.Fqy(res_coeffs)

    ## Representative ##
    def representative(self, omega):
        """
        Construct a representative phi of a type. A new level is added with
        phi and V.
        """

        s = len(self.levels)
        lvl_s = self.levels[s-1]
        ef = lvl_s.e * lvl_s.f
        u = ef * lvl_s.V

        if s+1 > 1:
            txp = lvl_s.inv_h * (u // lvl_s.e)
            twist = lvl_s.z^(txp)
        else:
            twist = self.levels[0].Fq(1)

        # Reductum from magma
        res_pol = twist * lvl_s.Fqy(lvl_s.res_pol.coefficients()[:-1])
        u += lvl_s.f * lvl_s.h
        phi0 = self.construct(s, res_pol, 0, u)
        print "%d. Repr. = %s + phi^ef" % (s, phi0,)

        phi0 = phi0 + lvl_s.phi^(ef)
        new_level = self.add_new_level(phi0, omega, u)

    ## Construct ##
    def construct(self, i, res_pol, s, u):
        """
        This routine constructs a polynomial phi0 with integer coefficients
        such that:
          - deg phi0 < m_i + 1 and y^nu*R_i(phi0)(y) = res_pol(y), where
            nu = ord_y(res_pol).
          - The non-negative integers s,u are the coordinates of teh left
            endpoint of a segment of slope -lvl_i.slope supporting
            N_i(phi0)
        """
        assert i <= len(self.levels), "i (%d) must be <= #levels (%d)" % (i, len(self.levels),)
        lvl_i = self.levels[i-1]
        assert res_pol.degree() < lvl_i.f, "res_pol is too large."
        assert u + s*lvl_i.slope >= lvl_i.f*(lvl_i.e*lvl_i.V + lvl_i.h), "the point (s, u) is too low"

        var = lvl_i.phi**lvl_i.e
        phi0 = 0
        height = u - res_pol.degree()*lvl_i.h
        if i == 1:
            for a in reversed(res_pol.coefficients()):
                lift = ZZ[x](list(a.polynomial()))

                phi0 = (phi0 * var) + (lift * self.prime^height)
                height = height + lvl_i.h
        else:
            raise NotImplemented, "Construct for i > 1 not implemented"

        phi0 = phi0 * lvl_i.phi^s

        return phi0




class MontesTypeLevel(object):

    def __init__(self, phi, omega, p):
        self.phi = phi
        self.prime = p

        self.e = None
        self.f = None
        self.h = None
        self.prod_e = 1
        self.prod_f = 1

        self.V = 0
        self.omega = omega
        self.cutting_slope = 0
        self.refinements = [ ]

        self.log_pi = vector([1, 0])
        self.log_phi = vector([0, 1])
        self.log_gamma = None

        self.Fq = None
        self.z = None
        self.Fqy = None

        self.slope = None
        self.res_pol = None

    def copy(self):
        copy = MontesTypeLevel(self.phi, self.omega, self.prime)

        copy.phi = self.phi
        copy.prime = self.prime

        copy.e = self.e
        copy.f = self.f
        copy.h = self.h
        copy.prod_e = self.prod_e
        copy.prod_f = self.prod_f

        copy.V = self.V
        copy.omega = self.omega
        copy.cutting_slope = self.cutting_slope
        copy.refinements = self.refinements

        copy.log_pi = self.log_pi
        copy.log_phi = self.log_phi
        copy.log_gamma = self.log_gamma

        copy.Fq = self.Fq
        copy.z = self.z
        copy.Fqy = self.Fqy

        copy.slope = self.slope
        copy.res_pol = self.res_pol

        return copy

    def __unicode__(self):
        phi = self.phi
        if self.slope is not None:
            slope = self.slope
        else:
            slope = '-'
        if self.res_pol is not None:
            res_pol = self.res_pol
        else:
            res_pol = '-'

        return '(%s, %s, %s)' % (str(phi), str(slope), str(res_pol))

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return self.__unicode__()

def montes(K, p, basis=False):

    if p.is_prime is False:
        raise Exception, "Error: p must be prime."
    f = K.defining_polynomial()
    if f.is_monic() is False:
        raise Exception, "Error: def.pol. of K must be monic."
    # TODO: Add requirements that coefficients of f are integers
    
    res_pol = PolynomialRing(FiniteField(p), 'y')(f)

    for factor in res_pol.factor():
        print "Analysing irred. factor modulo p: %s" % (str(factor[0]),)

        tt = MontesType(K, p, factor[0], factor[1])
        montes_main_loop(K, p, tt)

    
def montes_main_loop(K, p, tt):

    type_stack = [ tt ]
    while len(type_stack) > 0:
        tt = type_stack.pop()

        r = len(tt.levels)
        lvl_r = tt.rth_level()

        phiadic, quotients = phi_expansion(tt.pol,
                                           lvl_r.phi,
                                           lvl_r.omega)

        ## Phi-Newton Polygon
        sides, side_devs = tt.phi_newton_polygon(r, phiadic)

        length_N = lvl_r.omega
        index_N = -lvl_r.cutting_slope * (length_N*(length_N-1) // 2)
        if length_N == 1:
            print "Found a factor of depth r = %d\n", r-1
            raise NotImplemented, "No implementation of finilising types."

        # Create copies of tt for each side to use
        side_indexes_types = [(0, tt)]
        if len(sides) > 1:
            side_indexes_types += [(i, tt.copy()) for i in [1..len(sides)-1]]
        side_indexes_types.reverse()

        print side_indexes_types
        previous_h = 0
        for i, tt in side_indexes_types:
            #for i in reversed([0..len(sides)-1]):
            
            lvl_r = tt.rth_level()
            side = sides[i]
            print "Analysing side: %s" % (str(side),)

            lvl_r.h = - side[0].numerator()
            lvl_r.e = side[0].denominator()
            lvl_r.slope = - side[0]
            lvl_r.inv_h = inverse_mod(lvl_r.h, lvl_r.e)
            
            lprime = (lvl_r.inv_h*lvl_r.h - 1) // lvl_r.e
            new_pi = list( lvl_r.inv_h * lvl_r.log_phi - lprime*lvl_r.log_pi )
            new_pi.append(0)
            lvl_r.log_gamma = lvl_r.e*lvl_r.log_phi - lvl_r.h*lvl_r.log_pi

            E = ZZ(side[2][0] - side[1][0])
            H = ZZ(side[2][1] - side[1][1])
            index_N += E + previous_h + ((E*H-E-H+(E // lvl_r.e)) // 2)
            previous_h += H
            
            ## Residual polynomial
            res_pol = tt.residual_polynomial(r, side_devs[i])
            res_pol = res_pol.monic()
            factors = res_pol.factor()
            
            factors_types = [ (factors[0], tt) ]
            if len(factors) > 1:
                factors_types += [ (factors[i], tt.copy()) for i in [1..len(factors)-1] ]

            for factor, tt in factors_types:
                print "Analysing factor of the Res.Pol. %s" % (factor[0],)

                print tt
                lvl_r = tt.rth_level()

                omega = factor[1]
                lvl_r.res_pol = factor[0]
                lvl_r.f = lvl_r.res_pol.degree()

                ## Representative
                tt.representative(omega)

                ## TODO: Working here.


def phi_expansion(f, phi, omega):
    q = f
    coeffs = [ ]
    quos = [ ]
    for j in [0..omega]:
        q, r = q.quo_rem(phi)
        coeffs.append(r)
        quos.append(q)

    return coeffs, quos

