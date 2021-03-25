# CCA2
Computational Convex Analysis (CCA) numerical library for bivariate functions
Copyright (C) 2008-2021, Yves Lucet

This toolbox implements functions for convex operations on bivariate 
piecewise-linear cubic functions.

Here are the main files provided by the package. 
Filename | Description
-- | --
README.txt   | The present file
PLQVC.m      | The main class implementing convexity detection algorithm, plotting, data structure...
PLQVCTest.m  | Unit tests

## Unit Tests: in PLQVCTest.m

```MATLAB
>> runtests('PLQVCTest')
Running PLQVCTest
.......... .......... .......... .......... ......
Done PLQVCTest
__________


ans = 

  1×46 TestResult array with properties:

    Name
    Passed
    Failed
    Incomplete
    Duration
    Details

Totals:
   46 Passed, 0 Failed, 0 Incomplete.
   2.9043 seconds testing time.
```   

## Credits:

See copyright.txt for the list of contributors to this library.
