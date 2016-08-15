## TODO list

* Document regularizers. **Doing**
* Document solve routines.
* Implement tests on regularizers (also for matrix- and complex-variables). **Doing**
* Test algorithms on matrix- and complex-variable problems.
* Give some structure to algorithms output, such as

		x, output = solve(A, b, g, x0)
where `x` containts the solution Array, `output` contains other details about the solution process. (Just as an example, another idea is to do like other well established Julia solvers do).
* Similarly for the options? Instead of having 1000 arguments (and non uniform ones: the memory parameter of L-BFGS only affects the zerofpr solver, not ista nor fista).
* In the algorithms, the initial value of gamma should be chosen in a better way.
* Is `solve` a good name for the main routine? Maybe `fit`? Maybe `regls`?
* Also, give the option to avoid doing line-search on gamma at all.
* Link [Travis](https://travis-ci.org/) to repo (once public).
* Link [Coveralls](https://coveralls.io/) to repo (once public).

## Code optimization

* See if static typing helps and is used in the best way.
* See if we can speed up LBFGS two-loop-recursion implementing it in C.

## Regularizers to add

* Indicator of an affine subspace.
* Nuclear norm.
* Elastic net.
* Indicator of the probability simplex?
* Indicator of L1 norm ball
