# Numerical Exercises 4

The C++ code included in the `MDNVE_MCNVT` folder can be compiled and executed using the same procedure described in the main repository `README` file. To automate the equilibration process I also included three bash scripts:

- `run_eq.sh`: used to run only the equilibration procedure, results (temperature only) dumped into `./NonEq/`. Three flags allowed `-s`,`-l`,`-g`, each one specifies the requested phase (in order solid, liquid and gas).
- `run_sim.sh`: used to run the equilibration and simulation procedure, results dumped into `./Eq/`. Four flags allowed, the same three as the previous script and `-v`, which prints additional information.
- `clean.sh`: clean all the output files.

All the scripts can be run by using the command `sh script_name.sh -flags`.