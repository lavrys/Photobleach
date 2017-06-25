Photobleach Step Counting
===========

Simple Python Code for counting photobleach events and specifying their exact time.

Currently this repository holds the code as it existed when it was used for our Molecular Biology of the cell paper: http://www.molbiolcell.org/content/early/2016/09/19/mbc.E16-06-0404.full.pdf+html
plus some minor modifications to improve speed and usability.

This includes the core evaluator/estimator (i.e the cost function) and auxiliary routines, but not the iterator that optimizes parameters, which is part of the Machine Learning shell. 

A Jupyter notebok is available, containing a full presentation/description of the scientific problem (at the freshman undergraduate level) and explaining the use and function of this code. It is possible to run your data through the notebook, but many code parameters (preset in the notebook) could give you a better result if you changed their value, which at this time you can only do by editing the code and running it from the command line.

A full Machine Learning version with this code at its core is described in a publication which is currently under peer review. As soon as the review is complete the ML version will be published in this repository. This will happen irrespective of whether the review is positive or negative with regards to publication in the reviewing journal; since I am no longer in academia I do not really care that this be published in a journal, so if it is rejected I wil just put it here. 

Konstantinos Tsekouras

