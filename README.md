# Remote.Sensing
Hyperspectral image unmixing.

## create4.m getSynData.m
Code for generating spectral samples with dirichlet distribution.

## hyperVca.m
VCA implementation. testVCA.m is the corresponding test code.

## Nfindr.m
Nfindr implementation. testNfindr.m is the corresponding test code.

## hyperNmfMDC.m
NMF with minimum distance constraint on endmembers.

## hyperNmfMVC.m
NMF with minimum volume constraint on endmembers.

## hyperNmfASCL1.m
NMF with L1 sparsity constraint on abundance matrix.

## hyperNmfASCL1_2.m
NMF with L1_2 sparsity constraint on abundance matrix. L1_2 sparsity constraint seems inconsistent with reconstruction error. 

During the test, if constraint of sum-to-one constraint is set to be small (less than 1), L1_2 constraint will get very small abundance matrix (sum of coefficient at one pixel is far smaller than 1). However, if constraint of sum-to-one is set to me large (5~10 in synthetic data), resultant endmembers will be close to true value. 

Another thing is: when the value of reconstruction error and L1_2 constraint is in the similar level, reconstruction error will first decrease, while L1_2 constraint normally keep unchanged or slightly increase. when reconstruction error go down to a value far smaller than L1_2 constraint, L1_2 constraint will be the dominant factor in the following iterations. Namely, decrease of reconstruction error and L1_2 constraint are not in the same cycles. In some situations in which constraint of sum-to-one is weak, reconstruction error will go up again during decreasing of L1_2 constraint.

## test*
test codes for algorithms.

## exp*
detailed experiments.

