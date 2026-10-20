# iodine-xenon-135-dynamics

Python model of iodine-135 / xenon-135 buildup in a PWR core at constant power, solved with SciPy. The script also contains a standalone 6-group point kinetics solver for a step reactivity insertion.

## What it does

**Block 1: point kinetics (6 delayed neutron groups)**
Integrates the point kinetics equations for a 0.5 $ reactivity step, starting from equilibrium. The power plot is commented out in the script; uncomment it to display it.

**Block 2: I-135 / Xe-135 dynamics**
Integrates the two coupled balance equations, starting from zero concentrations:

```
dI/dt  = gamma_I * R - lambda_I * I
dXe/dt = gamma_Xe * R + lambda_I * I - lambda_Xe * Xe - sigma_a * phi * Xe
```

where `R` is the fission rate density (fissions·cm⁻³·s⁻¹) and `phi` the thermal neutron flux (n·cm⁻²·s⁻¹).

The two blocks are independent: the xenon block does not use the output of the kinetics block.

## Parameters (xenon block)

| Parameter | Value | Note |
|---|---|---|
| `gamma_Iodine` | 0.063 | Cumulative I-135 yield (includes Te-135 precursor), about 6 % |
| `gamma_Xenon` | 0.0024 | Independent Xe-135 yield, about 0.25 % |
| `lambda_Iodine` | 2.90e-5 s⁻¹ | I-135 half-life about 6.6 h |
| `lambda_Xenon` | 2.10e-5 s⁻¹ | Xe-135 half-life about 9.2 h |
| `sigma_a_Xenon` | 2.6e-18 cm² | Thermal absorption, about 2.6 million barns |
| `phi` | 3.0e13 n·cm⁻²·s⁻¹ | Representative PWR thermal flux |
| `R_fission` | 3.0e12 fissions·cm⁻³·s⁻¹ | About 100 kW/L, typical PWR core average |

`phi` and `R_fission` are representative orders of magnitude, not data from a specific plant.

## Sanity check

At equilibrium, both derivatives are zero:

```
I_eq  = gamma_I * R / lambda_I                                          = 6.5e15 atoms/cm³
Xe_eq = (gamma_Xe * R + lambda_I * I_eq) / (lambda_Xe + sigma_a * phi)  = 2.0e15 atoms/cm³
```

The numerical solution after 200 000 s (about 55.6 h) gives about 6.50e15 and 1.97e15 atoms/cm³, within 0.5 % of these values.

## Limitations

- Constant flux and constant power: the xenon concentration does not feed back on the flux.
- Starts from zero concentrations, so the plot shows the approach to equilibrium, not a shutdown.
- The post-shutdown xenon peak is not modeled yet (it can be obtained by restarting from equilibrium with `R_fission = 0` and `phi = 0`).
- Not validated against plant data.

## Sources

- Wikipedia, [Xenon-135](https://en.wikipedia.org/wiki/Xenon-135): half-lives of Xe-135 (9.14 h) and I-135 (6.58 h), absorption cross-section of Xe-135 (2 to 3 million barns; 2.65×10⁻¹⁸ cm² at a flux of 10¹⁴ n·cm⁻²·s⁻¹), order of magnitude of the I-135 and Xe-135 fission yields (about 6 %).
- *Measuring the fission yield of Xe-135 for applications*, [Springer](https://link.springer.com/article/10.1007/s10967-026-10794-3): ENDF/B-VIII.1 values for thermal fission of U-235, cumulative Xe-135 yield 6.54 %, independent yields 0.08 % (Xe-135) and 0.18 % (Xe-135m).
- Hawkings, Edwards, Olmstead (1971), [*Independent yields of mass-135 xenons in the thermal neutron fission of 233U, 235U, 239Pu, and 241Pu*](https://cdnsciencepub.com/doi/10.1139/p71-094), Can. J. Phys.: ratio of independent Xe-135 yield to cumulative I-135 chain yield, 0.0357 for U-235.
- [*Independent yield of 135I in the thermal neutron fission of 235U. The half-life of 135Te*](https://www.sciencedirect.com/science/article/abs/pii/0022190269900013): independent I-135 yield is 47 ± 2 % of the chain yield, Te-135 half-life 18 ± 2 s (justifies using the cumulative yield for `gamma_Iodine`).

`gamma_Iodine = 0.063` is derived, not read from a table: ENDF/B-VIII.1 cumulative Xe-135 yield (6.54 %) minus the independent Xe-135 and Xe-135m yields (about 0.26 %) gives about 6.3 %.



