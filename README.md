FML
========

<div align = 'left'>
  <img src="fml.jpg" width="200" />
</div>

This is a repository for a tutorial on computing a fully marginalized likelihood with importance sampling. Methodological and science details can be found at [this arXiv link](http://arxiv.org/abs/1504.07995).


Prerequisites
-------------

The tutorial uses common modules from a standard Anaconda installation. There are a few other modules you will have to install separately through e.g. pip. Those are

```bash
pip install corner
```

for creating corner plots,

```bash
pip install emcee
```

for performing a Markov chain Monte Carlo for a specified dataset and model, and

```bash
pip install rebound
```

for integrating the planetary system and thus computing the model evaluations.


Attribution
-----------

Please cite [this](http://adsabs.harvard.edu/abs/2016MNRAS.455.2484N) if you use this method in your research.
