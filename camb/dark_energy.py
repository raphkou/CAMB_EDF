from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool
from scipy.interpolate import CubicSpline


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int),
        ("is_df_model", c_bool, "using the Dark Fluid model"),
        ("omch2_eff", c_double),
        ("Omega_DE_eff", c_double),
        ("Omega_c_eff", c_double),
        ("xi", c_double)
    ]

    def validate_params(self):
        return True


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0. The sound speed can also be parameterized as a function of the scale factor
    through set_cs2_a_table.

    See :meth:`.model.CAMBparams.set_initial_power_function` for a convenience constructor function to
    set a general interpolated P(k) model from a python function.

    """
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyEqnOfState'

    _fields_ = [
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("cs2", c_double, "fluid rest-frame sound speed squared"),
        ("use_tabulated_w", c_bool, "using an interpolated tabulated w(a) rather than w, wa above"),
        ("use_tabulated_cs2_a", c_bool, "using an interpolated tabulated cs2(a) rather than cs2 above"),
        ("use_tabulated_cs2_k", c_bool, "using an interpolated tabulated cs2(k) rather than cs2 above"),
        ("use_tabulated_cs2_ktau", c_bool, "using an interpolated tabulated cs2(k*tau) rather than cs2 above"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]

    _methods_ = [
        ('SetWTable', [numpy_1d, numpy_1d, POINTER(c_int)]),
        ('SetDeltaTable', [numpy_1d, numpy_1d, POINTER(c_int)]),
        ('SetCs2Table_a', [numpy_1d, numpy_1d, POINTER(c_int)]),
        ('SetCs2Table_k', [numpy_1d, numpy_1d, POINTER(c_int)]),
        ('SetCs2Table_ktau', [numpy_1d, numpy_1d, POINTER(c_int)]),
        ('grho_de', [POINTER(c_double)], c_double),
        ('grho_cdm', [POINTER(c_double)], c_double),
        ('w_de', [POINTER(c_double)], c_double),
        ('w_de_only', [POINTER(c_double)], c_double),
        ('dw_da', [POINTER(c_double), POINTER(c_int)], c_double)
    ]

    def set_params(self, w=-1.0, wa=0, cs2=1.0,
                   is_df_model=False, omch2_eff=0, ombh2=0, H0=67, xi=1,
                   amp_delta=None, amp_cs2=None, pars=None):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed squared
        """
        
        self.w = w
        self.wa = wa
        self.cs2 = cs2
        self.validate_params()
        
        if (is_df_model == False and amp_delta is not None):
            folder = "/Users/kou/Documents/Professionnel/Sussex/CAMB/data/"
            modes = np.load(folder+"Modes_Moss_pos.npy")
            a_i = np.load(folder+"a_i_Moss_pos.npy")
            Omega_i = amp_delta@modes[0:len(amp_delta)]
            
            pars.Omega_i = Omega_i
            pars.a_i = a_i
            
            Omega_m = pars.omegam
            G = 6.67430e-11
            kappa = 8*np.pi*G
            c = 2.99792458e8
            sigma_boltz = 5.670374419e-8
            Mpc = 3.085677581e22
            Omega_g = kappa/c**2*4*sigma_boltz/c**3*pars.TCMB**4*Mpc**2/(3*(pars.H0*1e3)**2)*c**2
            Omega_neutrino = 7./8*(4./11)**(4./3)*Omega_g*pars.N_eff
            Omega_r = Omega_g+Omega_neutrino
            Omega_Lambda = 1-Omega_m-Omega_r
            
            a = np.logspace(-7,0,250)
            beta = 6
            Omega_lcdm_i = Omega_m/a_i**3+Omega_r/a_i**4+Omega_Lambda
            a_vec = np.expand_dims(a,1)
            Omega_Moss_i = Omega_i*Omega_lcdm_i*(2*a_i**beta/(a_vec**beta+a_i**beta))**(6/beta)
            Omega_DE = Omega_Lambda + np.sum(Omega_Moss_i, axis=1)
            
            ln_Omega_DE = np.log(Omega_DE)
            ln_a = np.log(a)
            dln_Omega_DE = ln_Omega_DE[1:]-ln_Omega_DE[0:-1]
            dln_a = ln_a[1:]-ln_a[0:-1]
            w_DE = -1.-1./3.*dln_Omega_DE/dln_a
            w_DE = np.append(w_DE, w_DE[-1])
            
            self.set_w_a_table(a,w_DE)
        
        if (is_df_model == True):
            self.is_df_model = True
            self.omch2_eff = omch2_eff
            self.Omega_c_eff = omch2_eff/(H0/100)**2
            self.Omega_DE_eff = 1-(omch2_eff+ombh2)/(H0/100)**2
            self.xi = xi
            
            folder = "/Users/kou/Documents/Professionnel/Sussex/CAMB/data/"
            eigenvectors_delta = np.load(folder+"eigenvectors_delta_50_xi_0.1.npy")
            eigenvectors_cs2 = np.load(folder+"eigenvectors_cs2_50_xi_0.1.npy")
            a_i = np.load(folder+"a_i.npy")
            
            a = np.expand_dims(np.logspace(-7,0,250),1)
            Omega_i = eigenvectors_delta[:,0:len(amp_delta)]@amp_delta
            beta = 6
            delta_i = Omega_i*(2*a_i**beta/(a**beta+a_i**beta))**(6/beta)*(a/a_i)**3
            delta = np.sum(delta_i, axis=1)

            log_cs2_i = eigenvectors_cs2[:,0:len(amp_cs2)]@amp_cs2
            log_cs2 = log_cs2_i*(2*a_i**beta/(a**beta+a_i**beta))**(6/beta)*(a/a_i)**3
            log_cs2 = np.sum(log_cs2, axis=1)
            cs2 = 10**log_cs2
            
            self.set_Delta_a_table(a, delta)
            self.set_cs2_a_table(a, cs2)
            
            if pars is not None:
                pars.DF_a = a
                pars.DF_delta = delta
                pars.DF_cs2 = cs2
    

    def validate_params(self):
        if not self.use_tabulated_w and self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cublic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w):
            raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy w(a) arrays must end at a=1')
        if np.any(a <= 0):
            raise ValueError('Dark energy w(a) table cannot be set for a<=0')

        a = np.ascontiguousarray(a, dtype=np.float64)
        w = np.ascontiguousarray(w, dtype=np.float64)

        self.f_SetWTable(a, w, byref(c_int(len(a))))

        return self
    
    def set_Delta_a_table(self, a, delta):
        """
        Set delta(a) from numerical values (used as cublic spline). Note this is quite slow.

        :param a: array of scale factors
        :param delta: array of delta(a)
        :return: self
        """
        if len(a) != len(delta):
            raise ValueError('Dark energy delta(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy delta(a) arrays must end at a=1')
        if np.any(a <= 0):
            raise ValueError('Dark energy delta(a) table cannot be set for a<=0')

        a = np.ascontiguousarray(a, dtype=np.float64)
        delta = np.ascontiguousarray(delta, dtype=np.float64)

        self.f_SetDeltaTable(a, delta, byref(c_int(len(a))))

        return self

    def set_cs2_a_table(self, a, cs2):
        """
        Set cs2(a) from numerical values (used as cublic spline). 

        :param a: array of scale factors
        :param cs2: array of cs2(a)
        :return: self
        """
        if len(a) != len(cs2):
            raise ValueError('Dark energy cs2(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy cs2(a) arrays must end at a=1')
        if np.any(a <= 0):
            raise ValueError('Dark energy cs2(a) table cannot be set for a<=0')

        a = np.ascontiguousarray(a, dtype=np.float64)
        cs2 = np.ascontiguousarray(cs2, dtype=np.float64)

        self.f_SetCs2Table_a(a, cs2, byref(c_int(len(a))))

        return self
    
    def set_cs2_k_table(self, k, cs2):
        """
        Set cs2(k) from numerical values (used as cublic spline). 

        :param k: array of wavenumbers
        :param cs2: array of cs2(k)
        :return: self
        """
        if len(k) != len(cs2):
            raise ValueError('Dark energy cs2(k) table non-equal sized arrays')

        k = np.ascontiguousarray(k, dtype=np.float64)
        cs2 = np.ascontiguousarray(cs2, dtype=np.float64)

        self.f_SetCs2Table_k(k, cs2, byref(c_int(len(k))))

        return self
    
    def set_cs2_ktau_table(self, ktau, cs2):
        """
        Set cs2(k) from numerical values (used as cublic spline). 

        :param k: array of wavenumbers
        :param cs2: array of cs2(k)
        :return: self
        """
        if len(ktau) != len(cs2):
            raise ValueError('Dark energy cs2(k*tau) table non-equal sized arrays')

        ktau = np.ascontiguousarray(ktau, dtype=np.float64)
        cs2 = np.ascontiguousarray(cs2, dtype=np.float64)

        self.f_SetCs2Table_ktau(ktau, cs2, byref(c_int(len(ktau))))

        return self

    def __getstate__(self):
        if self.use_tabulated_w:
            raise TypeError("Cannot save class with splines")
        return super().__getstate__()


def update_DF_model(pars, H0):
    pars.DarkEnergy.omch2_eff = pars.omch2_eff
    pars.DarkEnergy.Omega_c_eff = pars.omch2_eff/(H0/100)**2
    pars.DarkEnergy.Omega_DE_eff = 1-(pars.omch2_eff+pars.ombh2)/(H0/100)**2
    pars.DarkEnergy.set_Delta_a_table(pars.DF_a, pars.DF_delta)
    pars.H0 = H0
    
def update_DE_Moss_model(pars, H0):
    Omega_m = (pars.omch2+pars.ombh2)/(H0/100)**2
    G = 6.67430e-11
    kappa = 8*np.pi*G
    c = 2.99792458e8
    sigma_boltz = 5.670374419e-8
    Mpc = 3.085677581e22
    Omega_g = kappa/c**2*4*sigma_boltz/c**3*pars.TCMB**4*Mpc**2/(3*(H0*1e3)**2)*c**2
    Omega_neutrino = 7./8*(4./11)**(4./3)*grhog*pars.N_eff/(3*(H0*1e3)**2)*c**2
    Omega_r = Omega_g+Omega_neutrino
    Omega_Lambda = 1-Omega_m-Omega_r

    a = np.logspace(-7,0,250)
    beta = 6
    Omega_lcdm_i = Omega_m/a_i**3+Omega_r/a_i**4+Omega_Lambda
    a_vec = np.expand_dims(a,1)
    Omega_Moss_i = pars.Omega_i*Omega_lcdm_i*(2*a_i**beta/(a_vec**beta+a_i**beta))**(6/beta)
    Omega_DE = Omega_Lambda + np.sum(Omega_Moss_i, axis=1)

    ln_Omega_DE = np.log(Omega_DE)
    ln_a = np.log(a)
    dln_Omega_DE = ln_Omega_DE[1:]-ln_Omega_DE[0:-1]
    dln_a = ln_a[1:]-ln_a[0:-1]
    w_DE = -1.-1./3.*dln_Omega_DE/dln_a
    w_DE = np.append(w_DE, w_DE[-1])

    pars.DarkEnergy.set_w_a_table(a,w_DE)
    pars.H0 = H0


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed single fluid model
    (as for single-field quintessense).

    """

    _fortran_class_module_ = 'DarkEnergyFluid'
    _fortran_class_name_ = 'TDarkEnergyFluid'

    def validate_params(self):
        super().validate_params()
        if not self.use_tabulated_w:
            if self.wa and (self.w < -1 - 1e-6 or 1 + self.w + self.wa < - 1e-6):
                raise CAMBError('fluid dark energy model does not support w crossing -1')

    def set_w_a_table(self, a, w):
        # check w array has elements that are all the same sign or zero
        if np.sign(1 + np.max(w)) - np.sign(1 + np.min(w)) == 2:
            raise ValueError('fluid dark energy model does not support w crossing -1')
        super().set_w_a_table(a, w)

    def set_cs2_a_table(self, a, cs2):
        # check cs2 array has positive elements
        if np.any(cs2<0):
            raise ValueError('fluid dark energy model does not support cs2<0')
        super().set_cs2_a_table(a, cs2)
        
    def set_cs2_k_table(self, k, cs2):
        # check cs2 array has positive elements
        if np.any(cs2<0):
            raise ValueError('fluid dark energy model does not support cs2<0')
        super().set_cs2_k_table(k, cs2)
        
    def set_cs2_k_table(self, ktau, cs2):
        # check cs2 array has positive elements
        if np.any(cs2<0):
            raise ValueError('fluid dark energy model does not support cs2<0')
        super().set_cs2_ktau_table(ktau, cs2)


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementating the w, wa or splined w(a) parameterization in the PPF perturbation approximation
    (`arXiv:0808.3125 <https://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specifc (early) dark energy fluid model
    (`arXiv:1806.10608 <https://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """
    _fields_ = [
        ("w_n", c_double, "effective equation of state parameter"),
        ("fde_zc", c_double, "energy density fraction at z=zc"),
        ("zc", c_double, "decay transition redshift (not same as peak of energy density fraction)"),
        ("theta_i", c_double, "initial condition field value")]

    _fortran_class_name_ = 'TAxionEffectiveFluid'
    _fortran_class_module_ = 'DarkEnergyFluid'

    def set_params(self, w_n, fde_zc, zc, theta_i=None, pars=None):
        self.w_n = w_n
        self.fde_zc = fde_zc
        self.zc = zc
        if theta_i is not None:
            self.theta_i = theta_i


# base class for scalar field quintessence models
class Quintessence(DarkEnergyModel):
    r"""
    Abstract base class for single scalar field quintessence models.

    For each model the field value and derivative are stored and splined at sampled scale factor values.

    To implement a new model, need to define a new derived class in Fortran,
    defining Vofphi and setting up initial conditions and interpolation tables (see TEarlyQuintessence as example).

    """
    _fields_ = [
        ("DebugLevel", c_int),
        ("astart", c_double),
        ("integrate_tol", c_double),
        ("sampled_a", AllocatableArrayDouble),
        ("phi_a", AllocatableArrayDouble),
        ("phidot_a", AllocatableArrayDouble),
        ("__npoints_linear", c_int),
        ("__npoints_log", c_int),
        ("__dloga", c_double),
        ("__da", c_double),
        ("__log_astart", c_double),
        ("__max_a_log", c_double),
        ("__ddphi_a", AllocatableArrayDouble),
        ("__ddphidot_a", AllocatableArrayDouble),
        ("__state", f_pointer)
    ]
    _fortran_class_module_ = 'Quintessence'

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")


@fortran_class
class EarlyQuintessence(Quintessence):
    r"""
    Example early quintessence (axion-like, as arXiv:1908.06995) with potential

     V(\phi) = m^2f^2 (1 - cos(\phi/f))^n + \Lambda_{cosmological constant}

    """

    _fields_ = [
        ("n", c_double, "power index for potential"),
        ("f", c_double, r"f/Mpl (sqrt(8\piG)f); only used for initial search value when use_zc is True"),
        ("m", c_double, "mass parameter in reduced Planck mass units; "
                        "only used for initial search value when use_zc is True"),
        ("theta_i", c_double, "phi/f initial field value"),
        ("frac_lambda0", c_double, "fraction of dark energy in cosmological constant today (approximated as 1)"),
        ("use_zc", c_bool, "solve for f, m to get specific critical reshift zc and fde_zc"),
        ("zc", c_double, "reshift of peak fractional early dark energy density"),
        ("fde_zc", c_double, "fraction of early dark energy density to total at peak"),
        ("npoints", c_int, "number of points for background integration spacing"),
        ("min_steps_per_osc", c_int, "minimumum number of steps per background oscillation scale"),
        ("fde", AllocatableArrayDouble, "after initialized, the calculated background early dark energy "
                                        "fractions at sampled_a"),
        ("__ddfde", AllocatableArrayDouble)
    ]
    _fortran_class_name_ = 'TEarlyQuintessence'

    def set_params(self, n, f=0.05, m=5e-54, theta_i=0.0, use_zc=True, zc=None, fde_zc=None):
        self.n = n
        self.f = f
        self.m = m
        self.theta_i = theta_i
        self.use_zc = use_zc
        if use_zc:
            if zc is None or fde_zc is None:
                raise ValueError("must set zc and fde_zc if using 'use_zc'")
            self.zc = zc
            self.fde_zc = fde_zc


# short names for models that support w/wa
F2003Class._class_names.update({'fluid': DarkEnergyFluid, 'ppf': DarkEnergyPPF})
