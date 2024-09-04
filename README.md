Solving Linear DSGE Models with Bernoulli Methods

The directory structure is as follows:

    "algorithm" contains the algorithms with bernoulli_matrix_quadratic.m being the main algorithm

    "improvement_experiment" contains the Numerically Problematic Parameterization experiment in Table 2. Run mmb_loop_new.m to reproduce and then disp_results.m to display the results

    "mmb_experiment" contains the MMB experiment starting at the zero matrix in Table 3 and Smets and Wouters (2007) comparison in Table 1. Run mmb_loop_new_5.m to reproduce and then plot_results_new.m to display the results of the MMB experiment and SM07_results.m to display the Smets and Wouters (2007) comparison of Table 1

    "mmb_experiment_2" contains the MMB experiment starting at the QZ solution in Table 4. Run mmb_loop_new.m to reproduce and plot_results_new.m to display the results

    "mmb_replication" contains the 99 models of the MMB used in the MMB experiments

    "policy_experiment_2" runs the grid size experiment for Smets and Wouters (2007) of Figure 1. Run policy_loop_fill_2_new_errors.m to reproduce and create_plots_2_new_errors.m to plot the results.
