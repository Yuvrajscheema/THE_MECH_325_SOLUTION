import flat_belts as fb

knowns = {
    fb.d: 2,
    fb.D: 4,
    fb.C: 9*12,
    fb.n_d: 1750,
    fb.b: 6,
    fb.t: 0.05,
    fb.gamma: 0.035,
    fb.H_nom: 2,
    fb.Ks: 1.25,
    fb.nd: 1.0,
    fb.mu: 0.5,
    fb.g: 32.17,
}

results = fb.solve_belt(knowns)

print("Flat Belt Computed Values:")
for var, val in results.items():
    try:
        # Try converting to float for numeric display
        numeric_val = float(val)
        print(f"{var}: {numeric_val:.3f}")
    except (TypeError, ValueError):
        # If not numeric, just print the symbolic expression
        print(f"{var}: {val}")

