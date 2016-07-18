FML
========

This is a repository for a tutorial on computing a fully marginalized likelihood with importance sampling. Methodological and science details can be found at `this arXiv link <http://arxiv.org/abs/1504.07995>`_.


Prerequisites
-------------

The tutorial uses common modules from a standard Anaconda installation. There are a few other modules you will have to install separately through e.g. pip. Those are

    pip install corner

for creating corner plots,

    pip install emcee

for performing a Markov chain Monte Carlo to get a posterior sample of your model parameters given the data and some specified model, and
   
   pip install rebound

for integrating the planetary system and thus the model evaluations.


Attribution
-----------

Please cite `this <http://adsabs.harvard.edu/abs/2016MNRAS.455.2484N>`_ if you use this method in your research.

