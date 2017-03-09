# superinfection
This is code associated with the paper

Estimating the true prevalence of superinfection from limited datasets
by DB Reeves, AM Margaret, C Johnston, and JT Schiffer

The code is written to use the EM algorithm to infer the average richness of a viral strain in a population of trial participants. This is important because the presence of superinfection, or any individual being multiply infected by a different strain of the same virus implies that vaccines might be difficult to design. We show how to simulate populations in order to check the method in run_super, and also how to check what happens if the evenness parameter \alpha is misspecified in run_super_misspec_alpha. The main code is in superinfection.R and there are two other utilities.
