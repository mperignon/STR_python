from math import sqrt, log10, pow


def calculate_Reynolds_number(D, nu=1e-6, g=9.81, rho_w=1000., rho_s=2650.):
    
    R = (rho_s - rho_w) / rho_w
    Re = sqrt(R * g * D) * D / nu

    failMessage  = "This equation is only valid for Reynolds numbers "\
                    "below 2.6e+06. The calculated Reynolds number for "\
                    "this particle is %.2g" % Re
    assert Re <= 2.6e6, failMessage
    
    print Re
    
    return Re
        