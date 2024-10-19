import numpy as np
from scipy.optimize import root_scalar

def eos_core_15(p, T):
    """
    EOS for hcp iron - silicate (FeSi7) computes volume, density, thermal expansivity, heat capacity at constant P, 
    Grüneisen parameter, and isothermal incompressibility.
    """
    # Constants
    #----------------------------------------------------------------------------------------------------------
    T0 = 300.0  # Reference temperature [K]
    MFe = 55.845  # Molar mass [g/mol]
    MSi = 28.0855  # Molar mass [g/mol]
    V0 = 11.266*2  # A^3 per atom
    KT0 = 462.13 # Bulk modulus at reference [GPa]
    KTP = -0.23  # Pressure derivative of bulk modulus
    gamma0 = 1.9547  # Grüneisen parameter at reference
    gammaInf = 0.827  # Grüneisen parameter at infinite compression
    beta = 0.826  # Parameter for gamma
    theta0 = 44.574  # Einstein temperature at reference [K]
    a0 = 2.121e-4  # Anharmonic parameter
    m = 1.891  # Anharmonic exponent
    g = 1.339  # Parameter for thermal contribution
    R = 8.3144621  # Gas constant [J/mol/K]
    XFe = (93/55.845)/(93/55.845 + 7/28.0855)  # Molar fraction Fe-7Si
    XSi = 1.0 - XFe  # Molar fraction Fe-7Si
    Z = XFe*26.0 + XSi*14 # average Atomic number
    M0 = XFe*MFe + XSi*MSi  # Average molar mass [g/mol]

    # Derived quantities
    PFG0 = 0.10036e9 * (Z / V0) ** (5.0 / 3.0)
    c0 = -np.log(3.0 * KT0 / PFG0)
    c2 = 3.0 / 2.0 * (KTP - 3.0) - c0

    # Solve for zeta (V/V0)**(1/3)
    zeta = zetax(p, T, c0, KT0, c2, V0, T0, gamma0, gammaInf, beta, theta0, R, a0, m)

    # Volume and density
    zeta3 = zeta ** 3
    V = V0 * zeta3
    rho = M0 / V

    # Grüneisen parameter
    gamma = gammax(zeta3, gamma0, gammaInf, beta)

    # Einstein temperature
    theta = thetax(zeta3, theta0, gamma0, gammaInf, beta)

    # Heat capacity at constant volume
    Cv = Cvx(T, theta, R)

    # Isothermal incompressibility KT
    KT = np.exp(c0 - c0 * zeta) * KT0 * (5.0 + zeta * (-4.0 + 2.0 * c2 * (-2.0 + zeta) * (-1.0 + zeta) +
               c0 * (-1.0 + zeta) * (-1.0 + c2 * (-1.0 + zeta) * zeta))) / zeta ** 5

    # Thermal contributions to incompressibility
    KT += (-(gamma ** 2 * (T * Cvx(T, theta, R) - T0 * Cvx(T0, theta, R))) +
          (gamma * (1.0 - beta + gamma) + beta * gammaInf) * (Ethx(T, theta, R) - Ethx(T0, theta, R))) / V

    # Non-harmonic and electronic contributions
    KT -= (dEeaxdx(zeta3, T, R, m, a0) - dEeaxdx(zeta3, T0, R, m, a0)) / V0 + (Eeax(zeta3, T, R, m, a0) -
           Eeax(zeta3, T0, R, m, a0)) / V

    # Add non-harmonic contributions to Cv
    Cv += Cveax(zeta3, T, R, m, a0)

    # Thermal expansivity and Cp
    alpha = gamma * Cv / (KT * V)
    Cp = Cv * (1.0 + alpha * gamma * T) / (M0 * 1.0e-9)  # From J/mol/K to J/kg/K

    return np.array([1.e-3 * V, rho, alpha, Cp, gamma, KT])


def gammax(x, gamma0, gammaInf, beta):
    """ Grüneisen parameter """
    return (gamma0 - gammaInf) * x ** beta + gammaInf


def thetax(x, theta0, gamma0, gammaInf, beta):
    """ Einstein temperature """
    return np.exp(-(-1.0 + x ** beta) * (gamma0 - gammaInf) / beta) * x ** (-gammaInf) * theta0


def Ethx(T, theta, R):
    """ Harmonic contribution to thermal energy """
    return 3.0 * R * (theta / 2.0 + theta / (np.exp(theta / T) - 1.0))


def Cvx(T, theta, R):
    """ Heat capacity at constant volume for Einstein model """
    return 3.0 * R * theta ** 2 * np.exp(theta / T) / ((-1.0 + np.exp(theta / T)) ** 2 * T ** 2)


def Eeax(x, T, R, m, a0):
    """ Non-harmonic and electronic contribution to thermal energy """
    return 3.0 / 2.0 * R * m * a0 * x ** m * T ** 2


def Cveax(x, T, R, m, a0):
    """ Non-harmonic and electronic contribution to heat capacity """
    return 3.0 * R * m * a0 * x ** m * T


def dEeaxdx(x, T, R, m, a0):
    """ Derivative with respect to x of non-harmonic and electronic contribution to thermal energy """
    return 3.0 / 2.0 * R * m ** 2 * a0 * x ** (m - 1) * T ** 2


def zetax(p, T, c0, KT0, c2, V0, T0, gamma0, gammaInf, beta, theta0, R, a0, m):
    """ Solve for zeta = (V/V0)^(1/3) using root finding """

    def HolzapfelEinstein(zeta):
        zeta3 = zeta ** 3
        fun = (3.0 * np.exp(c0 * (1.0 - zeta)) * KT0 * (1.0 - zeta) *
               (1.0 + c2 * (1.0 - zeta) * zeta)) / zeta ** 5
        fun += Eeax(zeta3, T, R, m, a0) / (V0 * zeta3) - Eeax(zeta3, T0, R, m, a0) / (V0 * zeta3)
        fun += ((Ethx(T, thetax(zeta3, theta0, gamma0, gammaInf, beta), R) - Ethx(T0, thetax(zeta3, theta0, gamma0, gammaInf, beta), R)) *
                gammax(zeta3, gamma0, gammaInf, beta)) / (V0 * zeta3)
        return p - fun

    sol = root_scalar(HolzapfelEinstein, bracket=[0.5, 1.1], method='bisect')
    
    if not sol.converged:
        raise RuntimeError("Root finding did not converge for zeta")

    return sol.root

eos_core_15(300, 5000)