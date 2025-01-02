# A Stochastic Geometry Framework for RIS-assisted OFDM Cellular Networks

This package implements the stochastic geometry framework presented in the paper "A Stochastic Geometry Framework for Performance Analysis of RIS-assisted OFDM Cellular Networks" (https://arxiv.org/abs/2310.06754), providing tools for analyzing the performance of RIS-assisted cellular networks. Furthermore, nodes that serve secondary signal sources in addition to the primary signal source can be generalized beyond the application of RIS.

## Major building blocks
This package consists in three components, where the first one is building upon the other two:
* **Stochastic Geometry**: An analytical framework for assessing wireless network performance. It enables analytical expressions for network metrics, taking into accound the probabilistic distributions describng spatial randomness of nodes and propagation fading.
* **Laplace Transform**: A mathematical tool that captures the distributional properties of stochastic point processes. 
* **Numerical Integration**: Essential for evaluating complex analytical expressions, particularly those arising from the application of stochastic geometry and Laplace transforms.

### Misc functions
Additionally, basic functions not categorized above are placed in the Utils folder for the moment, with potential future refactoring:
* Parse JSON configurations and generate output files
* Pathloss functions model signal attenuation between nodes


## Examples 
To apply this framework, we present several examples abstracted from my project:
* Computing Coverage probability
* Computing Ergodic rate

## Documentation
We provide a comprehensive PDF documentation to facilitate long-term research and understanding.


# Research Applications
Beyond the example folder, we've integrated simulations related to my PhD project to enrich the codebase and promote open-source research. While the core functionality remains focused on providing open-source tools for stochastic geometry modeling in wireless networks, these additional simulations offer concrete examples of academic applications.

# Usage Guide
* Verify compatibility with GCC 12+ for optimal performance, particularly when utilizing parallel processing. 
* Enable OpenMP for accelerated computations by compiling with the (-fopenmp) flag.
* In numerical calculations, approximate infinite values with a suitably large finite number.
* Practical warning to prevent integration misuse (avoid singularity) and malfunctions (avoid diverging functions)