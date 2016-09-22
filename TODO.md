## TODO list

* Document solve routines.
* Similarly for the options? Instead of having 1000 arguments (and non uniform ones: the memory parameter of L-BFGS only affects the zerofpr solver, not ista nor fista).
* Is `solve` a good name for the main routine? Maybe `fit`? Maybe `regls`?
* Link [Travis](https://travis-ci.org/) to repo (once public).
* Link [Coveralls](https://coveralls.io/) to repo (once public).

## Code optimization

* Check those allocations resulting from the test, a significant amount of time is devoted to garbage collection.
* See if static typing helps and is used in the best way.
* See if we can speed up LBFGS two-loop-recursion implementing it in C.

