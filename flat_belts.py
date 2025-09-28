import sympy as sp

# --- Symbols ---
D, d, C = sp.symbols('D d C', positive=True)           # pulley diameters, center distance
n_d, n_D = sp.symbols('n_d n_D', positive=True)       # RPMs
V, VR = sp.symbols('V VR', positive=True)             # belt speed, velocity ratio
b, t, gamma = sp.symbols('b t gamma', positive=True)  # belt width, thickness, density
phi_d, phi_D, L = sp.symbols('phi_d phi_D L', positive=True)  # wrap angles, length
H_nom, H_d, Ks, nd, H = sp.symbols('H_nom H_d Ks nd H', positive=True) # horsepower
T, omega = sp.symbols('T omega', positive=True)       # torque, angular speed
w, W = sp.symbols('w W', positive=True)               # weight per foot
F1a, F2, Fi, Fc = sp.symbols('F1a F2 Fi Fc', positive=True)  # forces
mu, g = sp.symbols('mu g', positive=True)             # friction coefficient, gravity

# --- Equations ---
eq_VR = sp.Eq(VR, D/d)
eq_nD = sp.Eq(n_D, n_d * d/D)
eq_V = sp.Eq(V, sp.pi * d * n_d / 12)
eq_w = sp.Eq(w, 12 * gamma * b * t)
eq_W = sp.Eq(W, gamma * b * t * (12/1))
eq_phi_d = sp.Eq(phi_d, sp.pi - 2*sp.asin((D - d)/(2*C)))
eq_phi_D = sp.Eq(phi_D, sp.pi + 2*sp.asin((D - d)/(2*C)))
eq_L = sp.Eq(L, sp.sqrt(4*C**2 - (D - d)**2) + (D*phi_D + d*phi_d)/2)
eq_Hd = sp.Eq(H_d, H_nom * Ks * nd)
eq_H_link = sp.Eq(H, H_d)
eq_T = sp.Eq(T, 63025 * H_d / n_d)
eq_Fc = sp.Eq(Fc, W/g * (V/60)**2)
# Coupled forces
eq_force1 = sp.Eq(F1a - F2, 2*T/d)
eq_force2 = sp.Eq(Fi, (F1a + F2)/2 - Fc)

force_eqs = [eq_force1, eq_force2]

equations = [
    eq_VR, eq_nD, eq_V, eq_w, eq_W,
    eq_phi_d, eq_phi_D, eq_L, eq_Hd, eq_H_link,
    eq_T, eq_Fc
]

# --- Solver ---
def solve_belt(knowns: dict):
    """
    Solve flat-belt design equations for missing variables, including coupled forces.
    Returns numeric values.
    """
    # Copy knowns
    knowns = dict(knowns)
    
    # First, solve non-force equations iteratively
    unknowns = {s for eq in equations for s in eq.free_symbols if s not in knowns}
    solved = True
    while solved and unknowns:
        solved = False
        for eq in equations:
            eq_unknowns = [s for s in eq.free_symbols if s not in knowns]
            if eq_unknowns:
                sol = sp.solve(eq.subs(knowns), eq_unknowns, dict=True)
                if sol:
                    for k, v in sol[0].items():
                        knowns[k] = v
                        solved = True
        unknowns = {s for eq in equations for s in eq.free_symbols if s not in knowns}
    
    # Then, solve for coupled forces (F1a, F2, Fi) if not already given
    force_unknowns = {F1a, F2, Fi} - set(knowns.keys())
    if force_unknowns:
        # Substitute knowns into force equations
        subbed_eqs = [eq.subs(knowns) for eq in force_eqs]
        sol_forces = sp.solve(
            subbed_eqs,
            [F1a, F2, Fi],
            dict=True
        )
        if sol_forces:
            knowns.update(sol_forces[0])
    
    # Convert all values to floats
    results_numeric = {}
    for k, v in knowns.items():
        try:
            results_numeric[k] = float(v)
        except (TypeError, ValueError):
            results_numeric[k] = v
    # Compute belt dip if possible: dip = C^2 * W / (96 * Fi)
    if C in results_numeric and W in results_numeric and Fi in results_numeric:
        results_numeric['dip'] = results_numeric[C]**2 * results_numeric[W] / (96 * results_numeric[Fi])
    return results_numeric

# --- Exported symbols ---
__all__ = [
    "D", "d", "C", "n_d", "n_D", "V", "VR",
    "b", "t", "gamma",
    "phi_d", "phi_D", "L",
    "H_nom", "H_d", "H", "Ks", "nd",
    "T", "omega",
    "w", "W",
    "F1a", "F2", "Fi", "Fc",
    "mu", "g",
    "solve_belt"
]
