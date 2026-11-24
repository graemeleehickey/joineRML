# JMbdirect 0.1.1: installation failure

Both joineRML 0.4.8 and 0.4.9 runs failed while installing JMbdirect against FastJM 1.7.2:

```text
Error: object ‘survfitjmcs’ is not exported by 'namespace:FastJM'
Execution halted
ERROR: lazy loading failed for package ‘JMbdirect’
```

JMbdirect imports this function from FastJM in its `NAMESPACE`. The same installation failure occurs with both joineRML versions and prevents subsequent checks. This is an existing dependency incompatibility, not evidence of a joineRML 0.4.9 regression. JMbdirect's compatibility remains unverified until it can be installed with compatible dependencies.

The complete installation output is retained in the [CRAN-version logs](logs/cran.log) and [candidate-version logs](logs/candidate.log). No dependencies were patched or downgraded for this comparison.
