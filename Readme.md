# Framework

## Major building blocks
This package comprises three key components
* **Stochastic Geometry**: A powerful framework for analyzing wireless network performance. It enables analytical expressions for network metrics, considering probabilistic distributions like spatial randomness and fading.
* **Laplace Transform**: A mathematical tool that efficiently captures the distributional properties of stochastic processes. 
* **Numerical Integration**: Essential for evaluating complex analytical expressions, particularly those arising from the application of stochastic geometry and Laplace transforms.

### Misc functions
Additionally, basic functions not categorized above are placed in the Utils module for now, with potential future refactoring:
* Parse JSON configurations and generate output files
* Pathloss functions model signal attenuation


## Examples 
To illustrate the practical utility of these building blocks, we present concrete examples derived from my project:
* Average interference
* Coverage probability
* Ergodic rate

## Documentation
We provide a comprehensive PDF documentation to facilitate long-term research and understanding.


# Research Applications
Beyond the example folder, we've integrated simulations related to my PhD project to enrich the codebase and promote open-source research. While the core functionality remains focused on providing open-source tools for stochastic geometry modeling in wireless networks, these additional simulations offer concrete examples of academic applications.

# Usage Guide
* Install dependencies: Boost (for one-dimensional numerical integration) and GCC (for parallel processing)
* Enable parallel processing with OpenMP using GCC compiler flags
* Practical guidance for handling infinite values (approximated in numerical computations)
* Practical Tips to Prevent Integration Misuse and Malfunctions