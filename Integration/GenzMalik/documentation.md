# Introduction

The Genz-Malik algorithm is an adaptive numerical integration algorithm for approximating integrals over n-dimensional rectangular regions. It is based on the van Dooren-de Ridder algorithm, but it uses a different integration rule and error estimation technique. For a wide spectrum of integrands, the Genz-Malik algorithm consistently demonstrates superior computational efficiency compared to alternative algorithms when n is in the range 2-8, *but the algorithm is sometimes unreliable[^1]*.

The Genz-Malik algorithm works by recursively dividing the integration region into smaller subregions until the desired accuracy is achieved. The algorithm uses an integration rule to approximate the integral over each subregion, and it uses an error estimation technique to determine which subregions need to be further subdivided.

[^1]: We first encountered the unreliable behavior of the integration algorithm when dealing with the function of the form $exp(-xy^4)$, where the rate of change along the two axes was significantly different. 


## Utilization 
To use the C++ function, follow these steps:

**template \<class F>**\
**static value_type integrate(F function, Cube ND_region, double expected_error, double estimated_error, int max_f_eval, int real_f_eval);**

* **F** is a reference to a function. Its type will be automatically determined from the context.
* **ND_region** refers to the n-dimensional space within which the integration operation will be carried out.
*  **expected_error** sets a tolerance for the maximum relative error in the computation. However, if the computational resources are not fully utilized, the actual error may exceed this specified limit.
* **estimated_error** is determined by calculating the difference between the 5th order approximation and the 7th order approximation.. 
* **max_f_eval** sets a limit on the number of function evaluations, effectively controlling the computational resources allocated to the process.
* **real_f_eval** is set the number of function call excuted once the integration process is finished.


### Reference:
