# Numerical Exercises 7

The C++ code included in the `MDNVE_MCNVT` folder can be compiled and executed using the same procedure described in the main repository `README` file. To automate the equilibration and simulation processes I also included five bash scripts:

- `run_eq.sh`: same as Lab4.
- `run_sim.sh`: same as Lab4.
- `clean.sh`: same as Lab4.
- `run_sim_autocorr.sh`: used to compute and dump the data necessary for the autocorrelation analysis. Same flags as `run_sim.sh`.
- `run_sim_ar.sh`: used to execute a simulation for Argon gas. Same flags as `run_sim.sh`, with an additional `-m` flag to use Monte Carlo algorithm.

All the scripts can be run by using the command `sh script_name.sh -flags`.
