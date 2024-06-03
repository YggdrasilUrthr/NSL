# Numerical Exercises 6

The C++ code included in the `ISING_1D` folder can be compiled and executed using the same procedure described in the main repository `README` file. To automate the simulation at different temperatures I also included two bash scripts:

- `run_sequence.sh`: runs a sequence of simulation at different temperatures separated by a fixed step (see the file content for the specific values). It accepts three ordered arguments:
    1. `h` parameter, defaults to 0.
    2. sampling method, defaults to Metropolis. Set this to 0 to use Gibbs sampling.
    3. list of space-separated numbers. These are the temperature for which the equilibration data will be dumped.
- `clean.sh`: clean all the output files.

All the scripts can be run by using the command `sh script_name.sh arguments`.