# Modeling of 2019–20 coronavirus pandemic data


- Application of the model to a few countries data, see [model.ipynb](https://github.com/egorsobolev/covid19/blob/master/model.ipynb)
- Comparision of the pandemic develompent in a few countries, see [countries.ipynb](https://github.com/egorsobolev/covid19/blob/master/countries.ipynb)

## Model

We applied [Gompertz curve](https://en.wikipedia.org/wiki/Gompertz_function) to model COVID-19 pandemic data:

$$f\left(t\right) = a\mathrm{e}^{-b\mathrm{e}^{-ct}}$$

where
- $a$ is an asymptote, since $\lim_{t \rightarrow \inf} a\mathrm{e}^{-b\mathrm{e}^{-ct}} = a\mathrm{e}^{0}=a$
- $b$ sets the displacement along the x-axis (translates the graph to the left or right).
- $c$ sets the growth rate (y scaling)

The rate of growth:

$$r=\frac{f'(t)}{f(t)}=bc\mathrm{e}^{-ct}$$

