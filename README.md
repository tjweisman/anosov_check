# Checking the straight-and-spaced condition for Anosov representations

This repository provides an example implementation of an algorithm for checking if a representation of a surface group in SL(n, R) is Anosov, described by [Max Riestenberg](https://sites.google.com/view/max-riestenberg/) in [[R](https://arxiv.org/abs/2409.08015)].

This code has been tested on Linux machines, and the instructions below assume that you are working in a Unix-like environment. Contact [Teddy Weisman](https://public.websites.umich.edu/~tjwei/) at [tjwei@umich.edu](mailto:tjwei@umich.edu) if you need assistance.

## Prerequisites and setup
 
You will need a working installation of Python 3.

The main prequisite for running this code is the [geometry_tools](https://public.websites.umich.edu/~tjwei/geometry_tools/geometry_tools.html) Python package. To install geometry_tools:

1. Download geometry_tools. If you have git installed, you can just clone the repository from GitHub by running `git clone https://github.com/tjweisman/geometry_tools.git`. Otherwise, go to the [geometry_tools GitHub page](https://github.com/tjweisman/geometry_tools), and download and extract the repository using the "Code" link.

2. (Optional.) Create a new virtual environment using e.g. [venv](https://docs.python.org/3/library/venv.html). 

3. From a terminal, activate your new virtual environment (if you created one), and run `pip install .` from the directory where you downloaded the repository. If you do not have a working `pip` you can run `python setup.py install` instead, but **this is not recommended.**

Once you have installed geometry_tools, you should also download the files in the anosov_check repository (if you haven't already). Either run `git clone https://github.com/tjweisman/anosov_check.git` or download and extract everything via the "Code" link.

## Running the algorithm

If you have [Jupyter](https://jupyter.org/) installed, you can run the `anosov_check.ipynb` notebook. Alternatively you can run the script `anosov_check.py` by running `python anosov_check.py` from the directory where you downloaded the files.

Both the Jupyter notbook and the script perform the straight-and-spaced check on a fixed Hitchin representation of a surface group in SL(3,R) (in fact, it is a representation in the Fuchsian locus), by iterating through all geodesic words of length 10 in the group (with respect to a fixed presentation).

On my personal laptop, the check takes approximately 20 minutes to run.

## Caveats

**This code performs all of its computations using hardware floating-point arithmetic, which means precision is limited and the results of the check are not rigorous.** This is a consequence of the fact that the geometry_tools package is built on top of numpy; numpy can do fast numerical calculations in C, but most of its linear algebra routines only work when operating on arrays of 64-bit floats. 

At a later date I may re-implement the algorithm using tools which support better precision and rigorous numerical checks, but this is a long-term plan.

## (Some) implementation details

At its core, running the straight-and-spaced algorithm is just a matter of iterating through every word of some fixed length N in a surface group, and performing some angle/distance estimates for each word. Since the group in question has exponential growth, doing this naively is slow. The implementation of the algorithm in this repository takes heavy advantage of a *finite-state automaton* which accepts a unique geodesic word for each element of the surface group. The automaton itself was produced by the [kbmag](https://gap-packages.github.io/kbmag/) program [HT23].

Part of the reason the automaton gives a significant speedup is that the operations performed on each word w of length N can be (very roughly) described as follows:
1. Split up w into a prefix w+ and a postfix w-, each of length N/2.
2. Do some computations on w+ and w-.
3. Compare the results of the w+ computations with the result of the w- computations to get some angle/distance estimates.

Using the automaton, it is possible to avoid performing the same "prefix" and "postfix" operations many times, which makes the process more efficient. The general idea is to do the following:

*At each state S of the automaton*:
1. Perform all of the "prefix" computations coming from words of length N/2 corresponding to paths ending at S.
2. Perform all of the "postfix" computations coming from words of length N/2 corresponding to paths starting at S.
3. Compare the results of all of the prefixes and postfixes from 1. and 2.

The algorithm still involves some redundancy, but this is still much faster than the naive approach.

## Credits

Code by Teddy Weisman, based on an initial implementation by Max Riestenberg

## References

[HT23] Holt, D. and GAP Team, T., kbmag, Knuth-Bendix on Monoids and Automatic Groups, Version 1.5.11 (2023)<br />
(Refereed GAP package), https://gap-packages.github.io/kbmag. 

[R] Riestenberg, J. M., Certifying Anosov representations, (2024)<br />
arXiv preprint, https://arxiv.org/abs/2409.08015.
