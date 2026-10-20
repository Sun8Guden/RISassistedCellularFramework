# Research drivers

The simulations behind the published figures, kept close to how they were run. The library in `include/risgeo` grew out of this code; these drivers are research scripts rather than examples of the library's style.

* `ris_performance/`: performance of RIS-assisted cellular networks (interference composition, ergodic rate vs. BS/RIS density, coverage holes, wedge and square RIS deployments). It has its own helpers (`LaplaceTransformFunctions/`, `LaplaceTransformPointProcess/`, ...), written before the library, and uses the library only for numerical integration.
* `techno_economic/`: derivatives of the coverage probability and ergodic rate w.r.t. the BS and RIS densities, used to compare the marginal value of deploying BSs vs. RISs. It is built on the library.

Each `simulation_*` folder holds a driver `main_*.cpp`, a short `description.md` and the `output.json` it produced. Build all drivers from the repository root with `make paper` (add `OPENMP=1` for parallel integration) or with CMake and `-DRISGEO_BUILD_PAPER=ON`. A driver writes `output.json` to the current directory, so run it from its own folder:

```bash
cd paper/ris_performance/simulation_alpha
../../../build/paper/ris_performance/simulation_alpha/main_alpha
```

`simulation_alpha` reproduces its committed `output.json` to machine precision. The other drivers are left in the configuration of their last run, so check the commented lines in each `main_*.cpp` before regenerating a figure.
