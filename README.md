Code needed to run Venus Spin-Climate Evolution Model (Venus-SCEM) to reproduce the results and main figures from: "No recent long-lived habitable states on Venus due to spin-climate interaction"

explanation of directories/files:

./model
    run_Venus_evolution.py: main model script. Takes in inputs for sensitivity tests (from venus_run_input). Sorts solutions by final spin state.

    venus_run_input: input parameters for main model script. See run_Venus_evolution.py for explainations of parameters. 

    run_Venus_many.sh: submits batch jobs 

    concatenate_sorting.py: merges all batch outputs into single .pkl file. Automatically called after batch jobs successfully finish. Cleans up run by deleting separate batch directories after merging.

    ./input_data
        albedo_Yang_matrix.npy: Numpy array of albedo vs. insolation and rotation rate taken from Yang et al. 2013. Interpolated to give albedo due to cloud cover in the habitable state.

        Q_withIT_330.npy: tidal quality factor from ocean tides reported in Green et al. 2019. To be interpolated in sensitivity test including ocean tidal dissipation.

        sig_330_Green.npy: corresponding forcing frequency in Green et al. 2019

To run the model submit run_Venus_many.sh [casename] [num_batches] [runs_per_batch]
ex:  ./run_Venus_many.sh baseline 50 1000
^ this will create a new run called "baseline" with 50 batches of 1000 runs each = 50,000 total simulations. 
User will need to edit the file to run on their machine (currently configured for Harvard FASRC cluster). 


./baseline

    all_solutuons.pkl: merged output of baseline case, 50,000 solutions. Used to make Figures 2-4 in main text.

    concatenate_all.out: out text file with summary of final state distribution.
    
./figures

    constants_functions.ipynb: containts constants and functions useful for plotting

    Figure_2.ipynb: creates Figure 2 of main text

    Figure_3.ipynb: creates Figure 3 of main text

    Figure_4.ipynb: creates Figure 4 of main text

./sensitivity_tests
    ./nu_1e-4: core viscosity = 1e-4 m2 s-1 (for CMF)
    ./nu_10:  core viscosity = 10 m2 s-1 (for CMF)
    ./Q_100: Constant-Q with Q = 100
    ./Q_30:  Constant-Q with Q = 30
    ./ps_0.5bar: habitable atmosphere ps = 0.5 bar
    ./ps_4bar: habitable atmosphere ps = 4 bar
    ./ps_10bar: habitable atmosphere ps = 10 bar
    ./SRG_280: runaway greenhouse limit set to ASR = 280 Wm-2
    ./SRG_320: runawya greenhouse limit set to ASR = 320 Wm-2
    ./qsteam_500: thermal tide in steam state qo = 500 Pa
    ./qsteam_4000: thermal tide in steam state qo = 4000 Pa
    ./OceanTides_Green_330mIT: including ocean tides in habitable state from reported Q in Green et al. 2019
    ./OceanTides_altQ: including ocean tides in habitable state using reported dissipation in Green et al. 2019 to    compute Q
    ./tanh_albedo_16day: using tanh function for albedo v. rotation, critical rotation = 16 days
    ./tanh_albedo_32day: using tanh function for albedo v. rotation, critical rotation = 32 days
    ./tanh_albedo_48day: using tanh function for albedo v. rotation, critical rotation = 48 days

^ all above directories contain the all_solutions.pkl file with solutions and the input parameter file (venus_run_input)
some also contain a different run_Venus_evolution.py file with source modifications

    ./sup_figures.py: python script for making the duration of habitable state vs. rotation rate and last habitable period for each sensitivity tests:
        python3 sup_figures.py [direc] [S_RG]
        python3 sup_figures.py qsteam_500 300 
        ^ makes the figure for qsteam_500 with the baseline S_RG = 300 Wm-2



