# Observation-conditioned thermal transfer deadlines for superconducting cables

This repository reproduces the continuous-uncertainty observer experiment from Dongyue Kan’s manuscript *Observation-Conditioned Thermal Feasibility of Delayed Load Transfer for Superconducting Cables*. It estimates how long a backup transfer can wait when conductor temperature is hidden, coolant temperature is measured with bounded error, and the transfer itself has a delay.

The experiment uses 12 synthetic confirmation configurations with 12 queries each. Four representations use the same coolant measurements and continuous prior. The primary comparison evaluates one box, six state cells, 36 state cells, and 36 joint state and heat-capacity cells. A separate 48-query stress check perturbs the heating and conductance assumptions.

## Run the experiment

Use Python 3.12 and run the commands from the repository root.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python src/continuous_observer.py development --smoke
```

The smoke check runs one development query and three representations. The reference deadlines are −1, 62, and 70 seconds for `single_box`, `state_bank`, and `parameter_bank`, respectively. A deadline of −1 records an abstention. Results appear in `results/smoke/`.

```bash
python src/continuous_observer.py development 2
python src/continuous_observer.py confirmation 2
python src/continuous_observer.py confirmation 1
python src/equal_budget_audit.py
python src/model_stress.py
python src/summarize.py
```

The commands use frozen inputs in `data/` and write under `results/`. Run the timing portion without concurrent workloads if comparing wall-clock measurements.

| Path | Contents |
| --- | --- |
| `src/thermal_model.py` | Synthetic two-temperature equations, loss surrogate, and balanced two-bus load model |
| `src/continuous_observer.py` | Observation updates, comparison trajectories, holding regions, and delayed-transfer replay |
| `src/equal_budget_audit.py` | The 36-cell state-only comparison and the timing protocol |
| `src/model_stress.py` | The misspecified-physics stress experiment |
| `data/` | Frozen synthetic queries, measurement errors, and protocol |
| `reference-results/` | Original per-query results, summaries, and the archived observer figure |

## Reference results and scope

In the archived 144-query confirmation, 36 state cells increased the mean nonnegative waiting deadline by 10.125 seconds relative to six cells. Twenty-five queries improved, the median paired gain was zero, and approvals increased from 27 to 34. The 36-cell joint state and capacity partition gave 27 approvals. Abstentions are mapped to zero only for mean and gain summaries. The archived 1-second refinement preserved all approval decisions.

The synthetic thermal constants and two-bus equations originate in the study’s software. The observer uses ordinary RK4 comparison trajectories and a separate adaptive SciPy replay of the synthetic plant. The archived protocol uses the label `DOP853` for this replay. The implementation calls `solve_ivp` with its default RK45 solver. This package preserves that behavior. The numerical results apply to the declared synthetic model and parameters. Field operation requires its own data and engineering validation.
