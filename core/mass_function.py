from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as np
import time

times = []
mass_function = None


def mass_function_power_law(a, m, m_max, beta):
    """
    Mass function based on double broken power law
    :param a:
    :param m:
    :param m_max:
    :param beta:
    :return:
    """
    out = a * np.power(m, beta + 1) * np.exp(-m / m_max)
    out /= np.power(m_max, beta + 1) * np.exp(-1)
    out /= np.sum(out)
    return out


def step(x1, x2, e):
    """
    Converts mass points to a exponential law
    :param x1:
    :param x2:
    :param e:
    :return:
    """
    m = [x1, x2]
    p = np.power(m, e)
    p /= p[0]
    return m, p


def mass_function_hcss(mass):
    """
    Mass function of a HCSS (broken power-law)
    :param mass: Mass points (important for the normalization
    :return: The probability of the different masses
    """
    m1, p1 = step(0.0008, 0.008, -0.3)
    m2, p2 = step(0.008, 0.5, -1.3)
    p2 *= p1[1]
    m3, p3 = step(0.5, 1, -2.3)
    p3 *= p2[1]
    m = []
    m.extend(m1)
    m.extend(m2)
    m.extend(m3)
    p = []
    p.extend(p1)
    p.extend(p2)
    p.extend(p3)

    inter = interp1d(np.log(m), np.log(p), fill_value='extrapolate',
                     bounds_error=False)

    out = inter(np.log(mass))
    out = np.exp(out)
    out /= np.sum(out)
    out = integrate(inter, mass)
    return out


def integrate(inter, masses):
    """
    Integrates the interpolated mass function
    :param inter: The interpolated mass function
    :param masses: THe mass steps
    :return: The probability of the different masses
    """
    int_mass = np.linspace(masses[0], masses[-1], num=masses.shape[0]*100)
    out = inter(np.log(int_mass))
    out = np.exp(out)
    inter = interp1d(int_mass, out, fill_value=0, bounds_error=False)

    prob = []
    m1 = [masses[0]]
    dm = (masses[1:]+masses[:-1])/2
    m1.extend(dm)
    m2 = []
    m2.extend(dm)
    m2.append(masses[-1])

    for sm, em in zip(m1, m2):
        prob.append(quad(inter, sm, em)[0])
    prob = np.array(prob)
    prob[0] *= 2
    prob /= np.sum(prob)
    return prob


def get_mass_function(mass):
    """
    Returns the current mass function
    :param mass: The mass steps
    :return: The mass function
    """
    if mass_function is None:
        return mass_function_hcss(mass)
    else:
        return mass_function(mass)


def get_masses(l, masses, prob):
    """
    Calculates the the amount of masses for HCSS with a statistical MC approach
    :param l: Number of bounded stars
    :param masses: Mass steps
    :param prob:
    :returns: Calculated masses
    """
    mass_sel = np.array([], dtype=np.int32)
    ms = []
    s = masses.shape[0]
    rand = np.random.rand
#    randint = np.random.randint
#    mass_sel2 = randint(s, size=l_s)
#    mass_prob = prob
    while len(ms) < l:
        rands = rand(s)
        
        dec = rands - prob
#        p = np.where(dec < 0)[0]
#        mass_sel = np.append(mass_sel, mass_sel2[p])
        ms.extend(masses[dec < 0])
    mass_sel = np.array(ms)
    mass_sel = mass_sel[:l]
#    mass_sel = masses[mass_sel]
    return mass_sel


def mass_counter(l, masses_r, masses, prop):

    mass_sel = get_masses(l, masses_r, prop)
    mass_sel = np.append(mass_sel, masses)
    return np.unique(mass_sel, return_counts=True)[-1]-1


def get_mass_sample(l, masses, prob):
    m_l = len(masses)
    masses_r = np.repeat(masses, 5*m_l)
    prob = np.repeat(prob, 5*m_l)
    return mass_counter(l, masses_r, masses, prob)


def test():
    l = 100
    masses = np.linspace(0.1, 2, num=191)
    masses = np.append(masses, np.linspace(0.8, 1.89, 1091))
    masses = np.unique(masses)
    prob = get_mass_function(masses)
    m = []
    for i in range(1000):
        print(i)
        t0 = time.time()
        m = get_mass_sample(l, masses, prob)
        times.append(time.time()-t0)
    print(m)
    print(np.sum(m[:50]))
    print(np.sum(m)/l)

    
if __name__ == '__main__':
    test()
