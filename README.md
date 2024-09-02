# SingleTrial Decoding/Encopding RSA
this repository contains several scripts simulation single-trial RSA introduced in Kikumoto & Mayr, 2020, PNAS

# Simulation for expanded space solution for possible false positives
As discussed in Kikumoto, Sameshima & Mayr, 2022, PsychScience (in the supplementary), when the RSA model of interest (e.g., conjunction) is equivalent to the target label for decoding, 
this method could generate false positives. Such false positive results tend to occur when multiple lower substitute representations (e.g., stimulus, response, and rule) co-exist. 

One solution is to express RSA models in the "expanded" task design space. For example, if there are factor A (3 levels) and factor B (3 levels) and we are interested in the conjunction of AxB, the vanila single-trial RSA would set the target label to be 9-way classification (AxB). Instead, we can add another factor C (2 levels) that we don't care to expand the task space (so 18-way classification). This means the conjunction model of interest (and other models, too) will be modified and it is no longer an identity matrix. Because the conjunction is no longer equivalent to the target label, it prevents false positives. 

There are two scripts for this simulation:
Conjunction_Simulate_F_Expand(3rules).R applies this logic to Exp.1 of Kikumoto & Mayr, 2020, PNAS (using even-odd trials as factor C)
Conjunction_Simulate_F_Expand(4rules).R applies this logic to Exp.2 of Kikumoto & Mayr, 2020, PNAS (here multiple conjunctions could be tested)


