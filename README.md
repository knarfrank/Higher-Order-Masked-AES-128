# Higher-Order-Masked-AES-128

Implemention in C of the higher-order masking scheme proposed in [0] with CPRR method from [1]. 

The AES implementation uses the Common Shares method [3], and the Random Reduction method [2] to increase the performance of the implementation.


### SecEvalCombined
The most interesting function is the 'SecEvalCombined' function which implements the Common Shares and Random Reduction from [2] and [3] respectively. 



### References


[0] "Provably Secure Higher Order Masking of AES"

[1] "Higher-Order Side Channel Security and Mask Refreshing"

[2] "Further Improving Efficiency of Higher-Order Masking Schemes by Decreasing Randomness Complexity"

[3] "Faster Evaluation of SBoxes via Common Shares"

