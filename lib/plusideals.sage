# -*- coding: utf-8 -*-

import pprint
from itertools import product
pp = pprint.PrettyPrinter(indent=4)

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

    slopes = [ Side(slope, list2point(cloud[a]), list2point(cloud[b])) for a, b, slope in slopes ]
    return slopes

def list2point(l):
    return Point(l[0], l[1])

class Point(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __unicode__(self):
        return '(%d, %d)' % (self.x, self.y,)

    def __repr__(self):
        return self.__unicode__()

    def list(self):
        return [self.x, self.y]

class Side(object):

    def __init__(self, slope, p1, p2):
        self.slope = slope
        self.p1 = p1
        self.p2 = p2

    def __unicode__(self):
        return '[%s, %s, %s]' % (unicode(self.slope), unicode(self.p1), unicode(self.p2),)

    def __repr__(self):
        return self.__unicode__()

    def width(self):
        return self.p2.x - self.p1.x 
    
    def height(self):
        return self.p1.y - self.p2.y

class MontesType(object):

    def __init__(self, K, p, varphi, omega, initial=True):
        self.parent = K
        self.pol = K.defining_polynomial()
        self.prime = p
        self.varphi = varphi
        self.levels = [ ]
        self.sfl = [0, 0, 0, 0] # Single Factor Lifting
        self.phiadic = [self.pol.parent()(0) for i in range(0, 4)]

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

        new_level.V = lvl_s.e*u
        new_level.prod_e = lvl_s.prod_e * lvl_s.e
        new_level.prod_f = lvl_s.prod_f * lvl_s.f
        
        # FIXME: p^f, f = what??
        # We want Fq to be the extension made by attaching lvl_s.res_pol to
        # lvl_s.Fqy, but that isn't possible in Sage just yet.
        new_level.Fq = FiniteField(self.prime^new_level.prod_f, 'w'+str(s))
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y'+str(s))

        new_level.embedding = Hom(lvl_s.Fq, new_level.Fq).list()[0]
         
        if lvl_s.f > 1:
            lifted_res_pol = new_level.Fqy(\
                    [new_level.embedding(c) for c in lvl_s.res_pol])

            for z in new_level.Fq.list():
                if lifted_res_pol(z) == 0:
                    new_level.z = z
                    break
            if new_level.z is None:
                raise Exception, "Error: could not find class of y in {0} (res_pol = {1})".format(new_level.Fq, lvl_s.res_pol)
        else:
            new_level.z = -list(lvl_s.res_pol)[0]

        print "Fq: %s, Fpy: %s, z: %s" % (str(new_level.Fq), str(new_level.Fqy), str(new_level.z),)
        
        self.levels.append(new_level)
        return new_level

    def refinement(self, more_factors=False):
        r = len(self.levels) - 1
        if more_factors is True:
            self.lvl(r).refinements.append([self.lvl(r).phi,
                                            self.lvl(r).slope])
        self.lvl(r).cuttingslope = ZZ(self.lvl(r+1).slope)
        self.lvl(r).phi = self.lvl(r+1).phi
        self.lvl(r).omega = self.lvl(r+1).omega

        self.remove_last_level()

    def last_level(self, phiadic, side, side_dev, last_psi=True):
        #tt.add_last_level(phiadic, sides[0], side_devs)

        r = len(self.levels)
        lvl_r = self.lvl(r)
        lvl_r.e = 1
        if r > 1:
            # FIXME: What does nur stand for?
            nur = sum([self.lvl(j).slope/self.lvl(j).prod_e for j in [1..r-1]])
            self.sfl = floor((lvl_r.V/lvl_r.prod_e)-nur)

        if side.p1.x == 0:
            slope = -side.slope
            lvl_r.h = ZZ(slope)
            # FIXME: Why do we just take the first 2 elements?
            self.phiadic[0:2] = phiadic[0:2]
            if last_psi:
                res_pol = self.residual_polynomial(r, side_dev)
                print "last monic res.pol. = %s (%s)" % (res_pol, res_pol.monic())
                lvl_r.res_pol = res_pol.monic()
                lvl_r.log_gamma = lvl_r.log_phi - (lvl_r.h * lvl_r.log_pi)
        else:
            slope = +Infinity
        lvl_r.slope = slope
        
        return lvl_r

    def remove_last_level(self):
        self.levels.pop()

    def lvl(self, i):
        ''' 1 indexed access to array levels.'''
        assert i > 0, 'tt.lvl(i) is 1 indexed, so i must be > 0.'
        return self.levels[i-1]

    def __unicode__(self):
        return '(%s; %s)' % (self.varphi,
                             '; '.join([unicode(l) for l in self.levels]))

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
            print "N_%d : v_%d(%s phi_%d) = %s + %d" % (i, i, phiadic[k], i, str(val), n)
            if abs(val) != Infinity:
                cloud.append([k, val + n])
                all_devs.append(dev)
            n += self.lvl(i).V

        sides = lower_convex_hull(cloud)
        abscissas = [ p[0] for p in cloud ]    

        for side in sides:
            height = ZZ(side.p1.y) # [1][1]
            dev = [ ]

            # [1][0], [2][0]
            for j in range(side.p1.x, side.p2.x + side.slope.denominator(), side.slope.denominator()):
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
                height += side.slope.numerator()
                
            dev.append(side.p1.list()) # list(side[1])
            side_devs.append(dev)    

        if len(sides) == 0:
            raise NotImplementedError, "Case for a zero side Newton Polygon has not been implemented"


        return sides, side_devs

    ## Value of polynomial a(x) in ZZ[x] ##
    def value(self, i, a):
        assert i <= len(self.levels)

        val = +Infinity
        if a == 0:
            if i == 1:
                devs = 0
            else:
                devs = []
            return val, devs

        if i == 1:
            val = min([valuation(c, self.prime) for c in list(a)])
            devs = a
        else:
            devs = []
            lvl_im1 = self.lvl(i-1)
            step = lvl_im1.V + lvl_im1.slope
            min_height = 0
            V_height = 0
            quot = a
            k = 0
            last = 0

            while quot != 0 and min_height < val:
                quot, ak = quot.quo_rem(lvl_im1.phi)
                new_val, dev = self.value(i-1, ak)
                candidate = new_val + min_height
                if candidate <= val:
                    if candidate < val:
                        val = candidate
                        first_abscissa = k
                        first_ordinate = new_val + V_height
                        devs = [ dev ]
                    else:
                        for j in range(last+lvl_im1.e*2, k+1, lvl_im1.e):
                            if i-1 == 1:
                                devs.append(0)
                            else:
                                devs.append([])
                        devs.append(dev)
                    last = k
                min_height += step
                V_height += lvl_im1.V
                k += 1
            # FIXME: Why are we so certain that first_abscissa and
            #        first_ordinate will have been set?
            devs.append([first_abscissa, first_ordinate])
            val = ZZ(lvl_im1.e * val)

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

        lvl_i = self.lvl(i)
        height = side_devs[-1][1]
        #print "height, side devs:", height, side_devs
        res_coeffs = [ ]

        for dev in side_devs[:-1]:
            if (i == 1 and dev == 0) or (i > 1 and len(dev) == 0): 
                res_coeffs.append(lvl_i.Fq(0))
            elif i == 1:
                # coefficients are polynomials too.
                coeff = ZZ[x](ZZ(dev) // self.prime^(height))
                res_coeffs.append(coeff(lvl_i.z))

                j = side_devs.index(dev)
                #print "%d. order Res.Pol. the easy way c_%d = (%d / %d) = %d (z_%d = %d)" % (i, j, dev, self.prime**height, dev // self.prime^(height), i-1, lvl_i.z)
            else:
                lvl_im1 = self.lvl(i-1)
                twist_exp = (dev[-1][0] - lvl_im1.inv_h*height) // lvl_im1.e
                #print "%d. order Res.Pol. twist exp: %s" % (i, twist_exp,)

                # coefficients are polynomials too.
                coeff = self.residual_polynomial(i-1, dev)
                lifted_coeff = lvl_i.Fqy([lvl_i.embedding(c) for c in coeff])
                print "------"
                print "coeff:", coeff, "parent:", coeff.parent()
                print "z:", lvl_i.z, "parent:", lvl_i.z.parent()
                print "------"
                res_coeffs.append(lvl_i.z^(twist_exp) * lifted_coeff(lvl_i.z))

                j = side_devs.index(dev)
                #print "%d. order Res.Pol. c_%d = %s (z_%d = %s)" % (
                #        i, j, lvl_i.z^(twist_exp)*coeff(lvl_i.z), i, lvl_i.z,)

            height = height - lvl_i.h

        res_pol = lvl_i.Fqy(res_coeffs)
        print "%d. order Res.Pol = %s" % (i, res_pol,)

        return res_pol

    ## Representative ##
    def representative(self, omega):
        """
        Construct a representative phi of a type. A new level is added with
        phi and V.
        """

        s = len(self.levels)
        lvl_s = self.lvl(s)
        ef = lvl_s.e * lvl_s.f
        u = ef * lvl_s.V

        if s > 1:
            txp = -self.lvl(s-1).inv_h * (u // self.lvl(s-1).e)
            twist = lvl_s.z^(txp)
        else:
            twist = lvl_s.Fq(1)

        # Reductum from magma
        res_pol = twist * lvl_s.Fqy(list(lvl_s.res_pol)[:-1])
        print "%s = %s * Reductum(%s) = %s * %s" % (res_pol, twist, lvl_s.res_pol, twist, lvl_s.Fqy(list(lvl_s.res_pol)[:-1]))
        u += lvl_s.f * lvl_s.h
        phi0 = self.construct(s, res_pol, 0, u)
        print "%d. Repr. = %s + phi^ef" % (s, phi0,)

        phi0 = phi0 + lvl_s.phi^(ef)
        new_level = self.add_new_level(phi0, omega, u)
        return new_level

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
        lvl_i = self.lvl(i)
        assert res_pol.degree() < lvl_i.f, "res_pol is too large."
        assert u + s*lvl_i.slope >= lvl_i.f*(lvl_i.e*lvl_i.V + lvl_i.h), "the point (s, u) is too low"

        var = lvl_i.phi^(lvl_i.e)
        phi0 = 0
        height = u - res_pol.degree()*lvl_i.h
        print "res pol coeffs: %s (deg = %d)" % (list(res_pol), res_pol.degree())
        print "%d. Construct res.pol. coeffs: %s (%s)" % (i, unicode(list(reversed(list(res_pol)))), res_pol)
        if i == 1:
            for a in reversed(list(res_pol)):
                lift = ZZ[x](list(a.polynomial()))
                print "eltseq (i=1) %s" % (list(a.polynomial()),)

                phi0 = (phi0 * var) + (lift * self.prime^height)
                height = height + lvl_i.h
        else:
            step = (lvl_i.e * lvl_i.V) + lvl_i.h
            new_V = u - (res_pol.degree() * step) - (s * lvl_i.V)
            lvl_im1 = self.lvl(i-1)
            for a in reversed(list(res_pol)):
                # FIXME: What does pj stand for?
                pj = 0
                if a != 0:
                    txp, s_im1 = divmod(lvl_im1.inv_h*height, lvl_im1.e)
                    u_im1 = (new_V - (s_im1 * lvl_im1.h)) // lvl_im1.e
                    c = (a*lvl_i.z^txp)
                    # Doing the equivalent of (if Eltseq existed)
                    #         lvl_im1.Fqy(Eltseq(c, lvl_im1.Fq))
                    if c.parent().base_ring() == lvl_im1.Fq:
                        new_res_pol = lvl_im1.Fqy(list(c.polynomial()))
                    else:
                        # FIXME: This is really ugly, we're searching all
                        # combinations of "coefficients" in a polynomial in z
                        # in order to find the correct one, something better
                        # should be done here.
                        eltseq = None
                        for cs in product(*[list(lvl_im1.Fq) for j in range(0, lvl_im1.f)]):
                            el = sum([lvl_i.embedding(cs[k])*lvl_i.z^k for k in range(0, lvl_im1.f)])
                            if el == c:
                                eltseq = list(cs)
                        if eltseq is None:
                            raise Exception, "Could not perform Eltseq({0}, {1})".format(c, lvl_im1.Fq)
                        print "eltseq (i=%d) %s" % (i, list(a.polynomial()),)
                        new_res_pol = lvl_im1.Fqy(eltseq)
                        pj = self.construct(i-1, new_res_pol, s_im1, u_im1)
                phi0 = (phi0 * var) + pj
                new_V += step
                height += lvl_i.h

        phi0 = phi0 * lvl_i.phi^(s)
        print "%d. Construct pol = %s" % (i, str(phi0))

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
        self.inv_h = None

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
        self.embedding = None

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
        copy.inv_h = self.inv_h

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
        copy.embedding = self.embedding

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

    print "Running Montes for {0} in Z_{1}".format(f, p)

    reps_OM = [ ]
    trees = [ ]
    
    res_pol = PolynomialRing(FiniteField(p), 'y')(f)

    for factor in res_pol.factor():
        print "Analysing irred. factor modulo p: %s" % (str(factor[0]),)

        tt = MontesType(K, p, factor[0], factor[1])
        tree = montes_main_loop(K, p, tt)

        reps_OM += tree
        trees.append(tree)

    if len(reps_OM) == 1:
        reps_OM[0].rth_level().phi = f
        reps_OM[0].rth_level().slope = +Infinity

    print "OM Representatives:"
    print len(reps_OM)
    for tt in reps_OM:
        print "   ", tt

    
def montes_main_loop(K, p, tt):
    
    total_index = 0
    leaves = [ ]
    type_stack = [ tt ]
    while len(type_stack) > 0:
        print "STARTING LOOP:\n  type stack:", type_stack
        tt = type_stack.pop()

        r = len(tt.levels)
        lvl_r = tt.rth_level()

        print "\n\n\nAnalysing type of order: %d (%s)"  % (r, unicode(lvl_r),)
        print "--------------------------------------------------------------------------------\n"

        phiadic, quotients = phi_expansion(tt.pol,
                                           lvl_r.phi,
                                           lvl_r.omega)

        ## Phi-Newton Polygon
        sides, side_devs = tt.phi_newton_polygon(r, phiadic)
        print "Sides of Newton polygon: %s" % (sides,)
        
        length_N = lvl_r.omega
        index_N = -lvl_r.cutting_slope * (length_N*(length_N-1) // 2)
        starting_prod_f = lvl_r.prod_f
        print "index N: %d, length N: %d" % (index_N, length_N)
        if length_N == 1:
            tt.last_level(phiadic, sides[0], side_devs[0])
            print "  #### Found a factor of depth r = %d (%s):\n    %s\n" % (r-1, lvl_r.phi, tt)
            leaves.append(tt)
            sides = []

        # Create copies of tt for each side to use
        if len(sides) > 0:
            side_indexes_types = [(0, tt)]
            if len(sides) > 1:
                side_indexes_types += [(i, tt.copy()) for i in [1..len(sides)-1]]
                more_factors = True
            else:
                more_factors = False
            side_indexes_types.reverse()
            print side_indexes_types
        else:
            side_indexes_types = [ ]

        previous_h = 0
        for i, tt in side_indexes_types:
            #for i in reversed([0..len(sides)-1]):
            
            lvl_r = tt.rth_level()
            side = sides[i]
            print "Analysing side: %s" % (str(side),)

            lvl_r.h = - side.slope.numerator()
            lvl_r.e = side.slope.denominator()
            lvl_r.slope = - side.slope
            lvl_r.inv_h = inverse_mod(lvl_r.h, lvl_r.e)
            
            lprime = (lvl_r.inv_h*lvl_r.h - 1) // lvl_r.e
            new_pi = list( lvl_r.inv_h * lvl_r.log_phi - lprime*lvl_r.log_pi )
            new_pi.append(0)
            lvl_r.log_gamma = lvl_r.e*lvl_r.log_phi - lvl_r.h*lvl_r.log_pi

            #E = ZZ(side[2][0] - side[1][0]) # [2][0] - [1][0]
            #H = ZZ(side[1][1] - side[2][1]) # [1][1] - [2][1]
            E = ZZ(side.width())
            H = ZZ(side.height())
            index_N += (E * previous_h) + ((E*H-E-H+(E // lvl_r.e)) // 2)
            print "E: %d, H: %d, side: %s, index N: %d" % (E, H, str(side), index_N,)
            previous_h += H
            
            ## Residual polynomial
            res_pol = tt.residual_polynomial(r, side_devs[i])
            res_pol = res_pol.monic()
            factors = res_pol.factor()
            
            factors_types = [ (factors[0], tt) ]
            if len(factors) > 1:
                factors_types += [ (factors[i], tt.copy()) for i in [1..len(factors)-1] ]
                more_factors = True

            for factor, tt in factors_types:
                print "Analysing factor of the Res.Pol. %s" % (factor[0],)

                print 'old_tt:', tt
                lvl_r = tt.rth_level()

                omega = factor[1]
                lvl_r.res_pol = factor[0]
                lvl_r.f = lvl_r.res_pol.degree()

                ## Representative
                # The new level (lvl_rp1) is already part of the type.
                lvl_rp1 = tt.representative(omega)
                print 'new_tt:', tt

                if lvl_r.phi.degree() == lvl_rp1.phi.degree():
                    # Non-optimal, refining level r
                    print 'Refining, cutting slope: %s' % (str(lvl_rp1.slope),)
                    tt.refinement(more_factors=more_factors)
                else:
                    # Proceeding to higher order
                    print 'Proceeding to higher order'
                    
                    tt.log_pi = vector(new_pi)
                    tt.log_phi = -(lvl_rp1.V // lvl_r.e) * lvl_r.log_pi
                    tt.log_phi = vector(list(tt.log_phi) + [1])

                # Push the new or refined type onto the stack.
                type_stack.append(tt)
                print "after appending to stack: %s" % (str(type_stack),)

            ## End of `factors' for loop
        ## End of `sides' for loop
        
        total_index += starting_prod_f * index_N
        print "Added %d * %d to the index (--> %d)" % (starting_prod_f, index_N, total_index,)

    return leaves

def phi_expansion(f, phi, omega):
    q = f
    coeffs = [ ]
    quos = [ ]
    for j in [0..omega]:
        q, r = q.quo_rem(phi)
        coeffs.append(r)
        quos.append(q)

    return coeffs, quos

