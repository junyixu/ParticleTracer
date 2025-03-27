## Tokamak fields

$$\mathbf{B} = \boldsymbol{\nabla} \xi + \frac{r^2}{qR} \boldsymbol{\nabla}\theta,$$
where $q$ is safety factor.

Jacobian:
$$J = \boldsymbol{\nabla} \xi^1 \cdot \boldsymbol{\nabla} \xi^2 \times \boldsymbol{\nabla} \xi^3 = \frac{1}{\sqrt{g}}$$

Gradient:
$$\boldsymbol{\nabla} = \boldsymbol{\nabla} \xi^i \frac{\partial}{\partial \xi^i}$$

Divergence:
$$\boldsymbol{\nabla} \cdot \mathbf{A} = \frac{1}{\sqrt{g}} \frac{\partial}{\partial \xi^i} ( A^i \sqrt{g} )$$

Curl:  
If $\mathbf{B} = \boldsymbol{\nabla} \times \mathbf{A}$ then
$$B^i =  \frac{1}{\sqrt{g}} \epsilon^{ijk} \frac{\partial A_k}{\partial \xi^j}$$

$$\sin \zeta = -\frac{y}{R} \quad \cos \zeta = \frac{x}{R}$$

$$\begin{aligned}
  \hat{\zeta} = -\sin \zeta \hat{x} - \cos \zeta \hat{y}, \\
  \hat{\theta} = - \sin \theta \hat{R} + \cos \theta \hat{z}\\
  \hat{\mathbf{R}} = \cos \zeta \hat{x} - \sin \zeta \hat{y}.
\end{aligned}$$

From the above equations, it can be concluded that
$$\hat{\theta} = - \sin \theta \cos \zeta \hat{x} + \sin \theta \sin \zeta \hat{y} + \cos \theta \hat{z}.$$
