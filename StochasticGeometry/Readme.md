# Wireless Network Performance Analysis  
Stochastic geometry models can provide ergodic performance by averaging over both fading and spatial randomness.

## Basics
### Interference characterization
To capture both the spatial locations and independent propagation of interference,
the aggregated interference can be analyzed by applying the Laplace transform to the marked point process. One can easily obtain the average value, quantile, etc, using the basic knowledge of characteristic funcions. 

### Coverage probability
The coverage probability can be expressed as the Laplace transform of the aggregate interference, noise, and reflected signals, assuming the direct link is subject to Rayleigh fading.

The case of assuming the direct link as a Gamma distribution (Nakagami-m fading) can be directly derived as a simple variant.

### Ergodic rate
This metric can be directly obtained from coverage probability and Shannon capacity (and its variants).

## Advanced topics 

### Separation of the positive part
The coverage probability expression is valid only when the aggregated interference and noise are positive, without the reflected signals in the classical framework. However, when considering reflected signals, the aggregated term can be negative, implying the coverage can be ensured by reflection alone. To address the issue of negative values, it is necessary to separate the positive and negative components of the Laplace transform.