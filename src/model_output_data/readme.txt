Readme file for model_output_data.
Model output data for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
Authors: Bethan Gregory, Rodney Elliott, Justin  Deighan, Hannes Groeller, and Michael Chaffin

- The directories 'LSA' (for low solar activity models) and 'HSA' (for high solar activity models) contain the Monte Carlo output for the escape calculations for hydrogen produced by HCO+ dissociative recombination (table in Figure 2 and Table S1). There are eleven directories in each: one for each model run, labeled 1--10. The 11th directory ('example') is the one plotted in Figure 2, but is not included in the calculations for the average escape.

- In the LSA directory there is one further directory. The directory 'traced_positions_LSA' contains output from a separate LSA (1e5 test particles, HCO+ DR) run, where the particle positions were tracked. This is used to produce the plot of particle end positions in Figure 1.

- The configuration files (corona3d_2020.cfg and Hot_H.cfg), the executable (corona3d_2020), and the output files are given for each run.

- The output files include escape (console.out), density, and EDF data.

- Also included in 'model_output_data' are the escape probability data for low (escape_probabilities_LSA.csv) and high (escape_probabilities_HSA.csv) solar activity, which are also given in Table S2.
