
# Material for the tutorial paper *dynConfiR: The R package for sequential sampling models of decision confidence*

This repository contains the data and code used in the paper: *dynConfiR: The R package for sequential sampling models of decision confidence* (Hellmann, S., Zehetleitner, M., & Rausch, M., 2024, Preprint on PsyArXiv (<https://osf.io/preprints/psyarxiv/e354s>) [doi: 10.31234/osf.io/e354s](https://doi.org/10.31234/osf.io/e354s)).

## Structure:

-   The folder **application_Law** contains data and code for the Example section of the paper:

    -   the data by Law & Lee (Ng et al., 2021) (together with a reamde for the data)
    -   the file **Script_fitSeqSampConfModels_Law.R**, which contains the code for model fitting and prediction,
    -   the files **IC_analysis1.R** and **IC_analysis2.R**, which contain code for analysis and visualisation of the quantitative model comparison
    -   the folder **plotscripts** contains some R-scripts to produce the figures in the manuscript
    -   the folders **autosave_all_standard_models** and **autosave_fixed_noiseparameters** contain logging files and saved fitting results for individual participants and models (with subfolders for each fitted model)
    -   **parfits.RData** and **parfits_fixed_and_free.RData** contain the fitted parameters for the standard models and the standard models together with the restricted models with fixed noise parameters, respectively
    -   **fitsandpredictions_1.RData** and **fitsandpredictions_allfits_2.RData** contain the fitted parameters as well as predicted distributions for the standard models and the standard models together with the restricted models with fixed noise parameters, respectively

-   The folders **density_precision_dynaViTE** and **density_precision_RM** contain the code for the Precision analyis section in the manuscript, each containing:

    -   a R-script **Script_density_precision...** containing the main code for the computations and the visualization of the results,
    -   helper scripts **paramerer_prior.R** and **helper_compute_probs.R** (each containing identical code for both folders), used to simulate parameter sets and to compute the densities for different precisions
    -   several **.RData**-files containing the intermediate and final results of the analysis

-   The folder **recovery_analysis** with code used for the Recovery analysis section of the manuscript. It contains:

    -   The script **Script_par_rec.R**, which contains the main recovery analysis, from simulating artificial data and fitting the models to the simulated data, it sources several functions from
    -   the folder **helper_fct.R**, which contains R-files each defining functions:
        -   to collect and combine previous parameter fits from empirical studies (**read_and_collect_previous_fits.R**, the result is saved in **collected_fits_models.RData**),
        -   to sample random proportions of confidence samples from a Dirichlet distribution (**fun_sample_rating_props.R**, the information for this is saved in **collected_rating_proportions_and_alpha.RData**)
        -   to simulate artificial data using the sampled parameter sets and confidence proportions (**helper_simulate_artificialdata_fivesteps.R**)
        -   and fitting a specific model to given data and combine with known true parameters (**fitting_function_parrec_fivesteps.R**)
    -   the folder **prevfits** including **.RData**-files that contain estimated parameters from model fitting to different empirical data sets (different files) for different models (contained in the same file)
    -   the folder **saved_details**, containing **all_par_samples_list.RData** (with all parameter sets used as true parameters for the recovery) and the folder **fit_results** containing subfolders for each model, each containing the fitting results for each iteration of the recovery simulation
    -   the script **gather_plot_results.R**, which gathers all the results within the folder **saved_data/fit_results**, and visualizes the results of the parameter recovery study

-   dynConfiR-source package file (.tar.gz). This is an intermediate 0.0.4-version of the package with likelihood and fitting functions, used for the computations. A more recent version is available on [GitHub](https://github.com/SeHellmann/dynConfiR) and on [CRAN](https://cran.r-project.org/web/packages/dynConfiR/index.html). The newer versions on GitHub and CRAN may produce different results. The source package may be installed via

    ``` R
    install.packages("dynWEV_0.0.4.tar.gz", type = "source", dependencies=TRUE,repos="http://a.cran.mirror") 
    ```
    The more recent versions may be installed using `install.packages("dynConfiR")`


## Dependencies and compatibility:
The analysis use different other R packages. The dependencies for the `dynConfiR`-package should be handled, when installing the package. The scripts use the package `rstudioapi`, which should be installed with any RStudio version, to find the script path for setting the working directory. The package `parallel` is used for parallelization but should also be installed in any R distribution. 

The following additional packages are used:
``` R         
install.packages(c("BayesFactor", "tidyverse", "ggh4x", "ggpubr", "scales"))
```

### References

Ng, L. C. H., Law, F. H. F., Lam, A. M. W., Or, C. C.-F., & Lee, A. L. F. (2021). Metacognitive adaptation revealed in serial dependence of visual confidence judgments. Journal of Vision, 21 (9), 2487. https://doi.org/10.1167/jov.21.9.2487 

## Contact

For comments, remarks, and questions please contact me: [sebastian.hellmann\@tum.de](mailto:sebastian.hellmann@tum.de)
