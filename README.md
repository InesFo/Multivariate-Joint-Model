# Multivariate Joint Model for Longitudinal and Time-to-event data: simulating data and fitting models

Consider the situation where we have, for several people:

1) a repeateadly-measured outcome over time (Y variable)
2) a time-to-event outcome (F variable).

For instance, Y can be the biomarker value and F can be the time since diagnostic until surgery. Because Y and F can the associated (e.g., surgery occurs eralier for those with higher biomarker values), it is necessary to jointly model them: [Y,F].

This model (Diggle et al., 2008) assumes that $Y$, the longitudinal process follows a normal distribution

$$\symbf{Y} \sim MVN(\symbf{\mu_Y},\symbf{V(\theta)})$$

It is common that F variables are right-skewed, so assume it follows a log-Normal distribution. Thus, a transformation of the time also results in a normal distribution:
$$S=\log(F) \sim N(\mu_S,\eta^2)$$

Thus, the joint distribution $[Y,S]$ is a multivariate normal distribution:

$$
\begin{bmatrix}
\symbf{Y} \\
S
\end{bmatrix}
\sim MVN
\begin{pmatrix}
\begin{bmatrix}
\symbf{\mu_Y} \\
\symbf{\mu_S} \\
\end{bmatrix}
,
\begin{bmatrix}
\symbf{V(\theta)} & \symbf{g(\phi)} \\
\symbf{g^T(\phi)} & \eta^2 \\
\end{bmatrix}
\end{pmatrix}
$$

<p align="center">
  <img width="217" height="240" src="https://user-images.githubusercontent.com/62517130/179856615-49b328b9-7dab-4356-a9b4-5598c024ed77.gif">
</p>

---

### Now... notation apart, how can we fit a multivariate joint model in R?  
### Using the lme() function!

<p align="center">
  <img width="240" height="200" src="https://user-images.githubusercontent.com/62517130/179857216-d67937c0-a740-4aa9-811b-9dea90df2aba.gif">
</p>

First, we need to specify a linear mixed model for $Y$ and for $S$:

$$Y_{ij} = \beta_0 + \beta_1*time_{ij} + U_i + Z^Y_{ij}$$

$$S_{i} = \beta_2 + U_i + Z^S_i$$

The key to implement the multivariate model [Y,S] is to transform it into a univariate model. For that, the two outcome variables are stacked, creating a new and unique variable, named, for instance, $value$. Then, binary dummy variables $DY$ and $DS$ are created to identify which values belong to the $Y$ and $S$ outcomes, respectively:

<div align="center">
|id   | time   | value      | DY         | DS        | variable  |
|---- |------- |----------- |----------- |---------- |---------- |
|1    | 1      | $y_{11}$   | 1          | 0         |Y          |
|1    | 2      | $y_{12}$   | 1          | 0         |Y          |
|1    | 3      | $y_{13}$   | 1          | 0         |Y          |
|2    | 1      | $y_{21}$   | 1          | 0         |Y          |
|2    | 2      | $y_{22}$   | 1          | 0         |Y          |
|2    | 3      | $y_{23}$   | 1          | 0         |Y          |
|3    | 1      | $y_{31}$   | 1          | 0         |Y          |
|3    | 2      | $y_{32}$   | 1          | 0         |Y          |
|3    | 3      | $y_{33}$   | 1          | 0         |Y          |
|1    | NA     | $s_{1}$    | 0          | 1         |S          |
|2    | NA     | $s_{2}$    | 0          | 1         |S          |
|3    | NA     | $s_{3}$    | 0          | 1         |S          |
</div>
  
Then, we can fit a single model that combines the two variables $Y$ and $S$: 

$$ value_{ij}= DY(\beta_0 + \beta_1 time_{ij} + U_i + Z^Y_{ij}) + DY(\beta_2 + U_i + Z^S_i) $$

which, in the lme() sintax equals to:
```
lme(value ~ 0 + DY + DS + DY:time , data = mydata , random = ~ 1 | id , weights = varIdent(form = ~1 |variable))
```
---
### What now?
To evaluate this approach, I simulated data based on the model and fitted the model. I did this for 1000 simulated datasets, to obtain a distribution of the estimated parameters. The simulation procedure followed:
<p align="center">
  <img width="302" height="440" src="https://user-images.githubusercontent.com/62517130/179860186-894dc864-a64b-40d3-bf24-5d5f0102d2c8.png">
</p>

At the end, we found that except for extremelly low parameter values, the model behaves nicely:

<p align="center">
  <img width="460" height="612" src="https://user-images.githubusercontent.com/62517130/179861143-973a1f7f-640f-4deb-91b0-5f4aa0c4d864.png">
</p>

<p align="center">
  <img src="https://user-images.githubusercontent.com/62517130/179861627-9175e8e9-26b8-40de-8874-a5bce89eacab.gif">
</p>

---
### Contents:




