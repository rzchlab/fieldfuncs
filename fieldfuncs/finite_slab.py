import numpy as np

def By_surface(x, w, d, j):
    """Magnetic field directed perpendicular to the
    surface of a conducting slab at the slabs surface,
    uniform current:
    
    Field at this top surface along +y (up)
    
         -----------------------------  y = d/2
         |                           |  
         |             Slab          | y = 0
         |                           |
         ----------------------------- y = -d/2
      x = -w/2        x = 0        x = w/2
      
    ^+y
    |
    |
    ---->+x
    
    z out of monitor towards keyboard
    
    Current density is uniform and along +z

    The equation is from this paper:
    "Distribution of the magnetic field induced by a current 
    passing through slabs in the superconducting and normal states"
    D.D. Prokof'ev
    Technical Physics June 2006, Volume 51, Issue 6, pp 675–682
    doi:10.1134/S1063784206060016
    
    Arguments:
        x (float): coordinate where field is computed (meters)
        w (float): width of slab in x direction (meters).
        d (float): depth of slab in y direction (meters).
        j (float): current density in slab (Amps/meters**2)
        
    Returns:
        By (float): Magnetic field at x (Tesla).
    """
    A = w - 2*x
    B = w + 2*x
    C = 2*d
    mu0_over_4pi = 1e-7
    return mu0_over_4pi * j * (-A * np.arctan(C / A) + 
                               B * np.arctan(C / B) + 
                               C/2 * np.log((B**2 + C**2) / (A**2 + C**2)))

def By_arb(x, y, w, d, j):
    """Magnetic field directed perpendicular to the
    surface of a conducting slab at an arbitrary point
    in the slab:
    
    Field at this point (the 'o') along +y (up)
              | 
              v    
         -----------------------------  y = d/2
         |    o                      |  
         |             Slab          | y = 0
         |                           |
         ----------------------------- y = -d/2
      x = -w/2        x = 0        x = w/2
      
    ^+y
    |
    |
    ---->+x
    
    z out of monitor towards keyboard
    
    Current density is uniform and along +z

    The equation is from this paper:
    "Distribution of the magnetic field induced by a current 
    passing through slabs in the superconducting and normal states"
    D.D. Prokof'ev
    Technical Physics June 2006, Volume 51, Issue 6, pp 675–682
    doi:10.1134/S1063784206060016
    
    Arguments:
        x (float): x coordinate where field is computed (meters).
        y (float): y coordinate where the field is computed (meters).
        w (float): width of slab in x direction (meters).
        d (float): depth of slab in y direction (meters).
        j (float): current density in slab (Amps/meters**2)
        
    Returns:
        By (float): Magnetic field at x (Tesla).
    """
    A = w - 2*x
    B = w + 2*x
    C = 2*d
    D = d - 2*y
    E = d + 2*y
    mu0_over_4pi = 1e-7
    atan = np.arctan
    ln = np.log
    return -mu0_over_4pi * j * (
          -A * (atan(D/A) + atan(E/A))
        +  B * (atan(E/B) + atan(D/B))
        +  y * ln((B**2 + D**2) * (A**2 + E**2) / (B**2 + E**2) / (A**2 + D**2))
        + d/2* ln((A**2 + E**2) * (A**2 + D**2) / (B**2 + D**2) / (B**2 + E**2))
    )
    

def By_2d_approximation(x, w, d, j):
    mu0_over_4pi = 1e-7
    return 2e-7 * j * d * np.log((w/2 + x) / (w/2 - x))

def g1d(x, x0, beam_fwhm):
    """1d Gaussian centered at x0 and with FWHM beam_fwhm."""
    s = beam_fwhm / 2.3548
    return twopi**-0.5 / s * np.exp(-(x - x0)**2/(2*s**2))


def box1d(x, x0, d):
    """Box function. Equal 1.0 for x less than d from x0 and 0 otherwise."""
    return (abs(x - x0) <=d).astype(int)


def get_f_amp(fwhm, w):
    """Create the parameterized beam smoothed amplitude function.
    
    Args:
        fwhm: Full Width Half Max of beam (in meters)
        w: width of metal strip (in meters)
    
    Returns:
        A function that can be used with scipy.optimize.curve_fit to
        fit data. The functions signature is f(x, center, height).
    """
    def f(x, center, height):
        xrange = 2 * fwhm
        y = []
        for xi in x:
            xx = np.linspace(xi - xrange, xi + xrange, 200)
            box = height * box1d(xx, center, w/2)
            y.append(simps(box * g1d(xx, xi, fwhm), xx))
        return np.array(y)
    return f

def get_f_off(fwhm, w, d, I):
    """Create the parameterized beam smoothed offset function.
    
    Args:
        fwhm: Full Width Half Max of beam (in meters)
        w: width of metal strip (in meters)
        d: thickness of metal strip (in meters)
        I: current passing through sample (in z direction)
    """
    j = I / (w * d)
    def f(x, center, height):
        xrange = 2 * fwhm
        y = []
        for xi in x:
            xx = np.linspace(xi - xrange, xi + xrange, 500)
            box = box1d(xx, center, w/2)
            By = height * By_surface(xx - center, w, d, j)
            y.append(simps(box * g1d(xx, xi, fwhm) * By, xx))
        return np.array(y)
    return f
