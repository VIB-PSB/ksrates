.. _`paralogs_analyses`:

Mixture modeling of paralog *K*:sub:`S` distributions
*****************************************************

The interpretation of mixed paralog–ortholog *K*:sub:`S` distributions is sometimes challenged by the fact that paralog WGD peaks are often not clearly distinguishable due to progressive WGD signal erosion over time and due to potential overlaps between peaks of successive WGDs. In order to more objectively define the *K*:sub:`S` age of WGD peaks, a clustering feature based on mixture modeling has been implemented in *ksrates*.

Three methods are available: anchor *K*:sub:`S` clustering, exponential-lognormal mixture modeling and lognormal-only mixture modeling.
For an extended description of these methods, please refer to the `Supplementary Materials <https://www.biorxiv.org/content/10.1101/2021.02.28.433234v1.supplementary-material>`__, currently in preprint. A short overview can be found below.

The analysis configuration determines which are applied and which are not (see also sections :ref:`pipeline_config_section` and :ref:`expert_config_section`).

    * If collinearity analysis is selected in the *ksrates* configuration file (``collinearity`` =  yes), the default method performed is a clustering based on the anchor pair *K*:sub:`S` values in collinear segment pairs. 
    * Otherwise (``collinearity`` = no and ``paranome`` = yes), the default method is exponential-lognormal mixture modeling of the whole-paranome *K*:sub:`S` distribution.
    * Lognormal-only mixture modeling is never applied by default, since it is more prone to produce spurious peaks.

The execution of non-default methods can be triggered in the expert configuration file and follows the table below. For example, when both collinearity and paranome analyses are selected (last column), only the anchor *K*:sub:`S` clustering is perfomed by default and all the other methods can be optionally triggered.

.. table:: Default methods are marked by a bold capitalized "X", while optional methods are marked by a lower-case "x".

    =======================================  ================  =============  ========================
    Method                                   Colinearity-only  Paranome-only  Collinearity and paranome
    =======================================  ================  =============  ========================
    Anchor *K*:sub:`S` clustering            **X**                            **X**
    Exponential-lognormal mixture model                        **X**          x
    Lognormal mixture model on anchor pairs  x                                x
    Lognormal mixture model on paranome                        x              x
    =======================================  ================  =============  ========================


.. _`anchor_ks_clustering`:

Anchor *K*:sub:`S` clustering
+++++++++++++++++++++++++++++

The anchor *K*:sub:`S` clustering is based on i-ADHoRe output files and on the ``.ks_anchors.tsv`` file produced by the pipeline. i-ADHoRe runs a colinearity analysis on the input sequences by using the GFF file, which contains the location of the sequences in the genome. When i-ADHoRe detects different regions that contain the same paralog genes in the same order (colinear blocks), it stacks and aligns such segments together to form a multiplicon. The genes in the multiplicon are called anchors and are likely to be paralogs due to large- or whole-genome multiplications. By performing a profile search, the program can build multiplicons that are composed of more than two segments, reflecting the presence of multiple WGMs. For more details refer to i-ADHoRe documentation. The *K*:sub:`S` associated to the anchor pairs in a multiplicon are called anchor *K*:sub:`S`, and they are estimated by ``wgd`` software and collected in the ``.ks_anchors.tsv`` file.

During the clustering, the attention is not placed on the multiplicon, but instead on all the possible segment pairs that the multiplicon is composed of. A segment pair originated from a recent WGD event will have low anchor *K*:sub:`S` values on it, while a segment pair originated from a more ancient WGD will show higher anchor *K*:sub:`S` values on it. The idea is to cluster the segment pairs found all over the multiplicons on the basis of their age.

Before clustering the segments, there is a redundancy problem that has to be taken into account. It often happens that two or more multiplicons share some anchor pairs. This means that there is redundancy among *K*:sub:`S` among all segment pairs, and that some anchor pairs/*K*:sub:`S` may be present in multiple copies. The filtering step compares all segment pairs and:

* segment pairs whose anchor pair list is a subset of the anchor pair list of another segment pairs are removed because absolutely redundant;

* if a segment pair has an anchor pair list that partially overlaps (> 1/3) with the list of another segment pair, the segment pair with the shortest list is removed.

* smaller overlapping is ignored.

Note that the filtering step considers the segment pair as a unit and not the single anchor pair in it because the colinearity traces are structured in blocks.

The age of a segment pair is ideally represented by the *K*:sub:`S` values of the anchor pairs contained in it: all its anchor *K*:sub:`S` values are expected to be similar because they were all generated at the same whole-genome duplications. However, there can be some variability in the *K*:sub:`S` estimates, with some values being outliers. If the segment pair contains few anchor pairs (up to 5), then its *K*:sub:`S` list is considered too short for reliable outliers detection; otherwise, then the *K*:sub:`S` falling outside the boundaries of the median absolute deviation are considered outliers and removed.

After that, for each segment pair it is computed the median of its anchor pair *K*:sub:`S` values to obtain a representative value of the segment pair age, which will be used in the clustering.

The clustering technique is Gaussian mixuture model from ``sklearn`` package and it first requires to set the number of clusters in advance. The number of segments in a multiplicon (its "level") provides a hint about how many WGMs happened in the species genome. However, it is not possible to know from the level itself if they were duplications or triplications, and in which combination; moreover, genome rearrangements and deletions can lead to a level that underestimates the number of multiplications. Our workflow simplifies this choice by assuming only duplications and no genome lost. It considers the highest level found among all multiplicons and infers the necessary minimum number of WGDs to explain it. For example, having the highest multiplicon level of 8 is explained by assuming 3 WGDs. The number of WGDs is given to the clustering function as the number of components.

A lognormal mixture model is then applied to the distribution of the segment pair *medians*. Then all *K*:sub:`S` lists belonging to medians clustered together are recalled, so to have clusters composed of segment pair *K*:sub:`S` lists and not anymore segment pair *K*:sub:`S` medians. This step produces a "raw" *K*:sub:`S` cluster plot.

Subsequently, a filtering step detects clusters with non-significant signal for peak calling and eliminates their *K*:sub:`S` values from the plot. One criterion among the following is enough for elimination: being poorly populated (containing 10% of total *K*:sub:`S` values), being very old (cluster median *K*:sub:`S` is >= 3 *K*:sub:`S`) and having quite a flat peak (inter-quartile range > 1.1 *K*:sub:`S`). The thresholds are based on *K*:sub:`S` plots of common species such as Arabidopsis, oil palm, rice.

Finally, a second round of lognormal mixture model is performed on the remaining *K*:sub:`S` values, using the number of clusters remaining after the filtering. Just as before, what is clustered are the medians of the remaining segment pairs. The clusters of the *K*:sub:`S` values are then again obtained from the clusters of medians. This second clustering step may look (slightly) different from the first clustering result since it starts from a different input median list.

The plot also includes the rate-adjusted divergence lines so to have a mixed plot where it is possible to compare the temporal relationship between the called WGD peaks and the speciation events. The final plot is saved in PDF format as ``mixed_species_anchor_clusters_unfiltered.pdf``, where species is the name of the focal species.


.. _`elmm`:

Exponential-lognormal mixture model
+++++++++++++++++++++++++++++++++++

A customized algorithm for mixture modeling with exponential and lognormal components has been implemented in the attempt to identify WGM traces in paranome distributions. The EM algorithm is based on Zhang et al. 2019. and is the default mixture modeling choice.

The exponential component is used to model the L-shaped background distribution of small-scale duplications (SSDs), whose chance to be kept tends to exponentially decrease as older they get. The lognormal components are instead used to model the WGM traces and are preferred to normal components because *K*:sub:`S` values can't be negatives and because WGMs tend to have a longer right tail.

The code performs the expectation-maximization (EM) algorithm to fit the mixture model on the paranome. Since the initialization of the component parameters plays an delicate role in mixture models, three strategies are followed and the best result is separately plotted: 1) guessing the parameters from the *K*:sub:`S` data itself, 2) starting with random parameters and 3) a hybrid initialization. In all the three strategies, an extra "buffer" lognormal component is by default initialized around 5 *K*:sub:`S` to avoid that the other components are forced to stretch towards higher values in the attempt to cover the entire distribution.

Despite the frequent application of this technique on *K*:sub:`S` distributions, the interpretation of mixture models requires caution due to the tendency of overfitting and overclustering. Note also that the reliability of components covering medium-high *K*:sub:`S` values (i.e. above 3 *K*:sub:`S`) is uncertain.


Initialization through data
---------------------------

The initialization of the rate of the exponential component is based on the height of the first bin in the paranome histogram plotted with bin width equal to 0.1 and the ``density`` option as ``True`` (namely, the area under the histogram integrates to 1, for more details see ``matplotlib`` documentation).

The initialization of the mean and standard deviation of the lognormal component(s) is based on the mean and standard deviation of the related normal curve(s). In order to guess this latter, the *K*:sub:`S` paranome is log-transformed and the new distribution is searched for peaks, which are potentially generated by putative WGM lognormal traces.

In more details: first, the KDE is obtained from the log-transformed data, then a spline is obtained from the KDE in order to smooth out the small noise irregularities on the KDE curve. The function find_peaks() from scipy package performs the peak detection on the spline. However, it can happen that some noise is mistaken for a peak or that some real peak signals are assigned a short prominence due to the overlapping with another close signal. To filter away the false positives and retain the real peaks, the distribution is mirrored in both directions around each peak and if the prominence of the peak is above an arbitrary threshold the peak is considered significant and retained. From the new prominences obtained after reflection it is guessed the width of the peak again through scipy. The peak x-coordinate is finally used as mean of the normal component and the peak width is used as standard deviation. If too wide, the standard deviation is reduced to a default intermediate value.

Random initialization
---------------------

The pure random choice of values can be quite misleading for the fitting or can require a large amount of initializations in order to obtain a good result. Therefore, the EM is initialized with components whose parameters are randomly taken from an appropriate arbitrary range of values for *K*:sub:`S` distributions. The fitting is performed considering the *K*:sub:`S` paranome histogram with ``density`` set to ``True`` (as in the other initialization method).

* The exponential rate is randomly chosen from a range between 0.2 and 1 with step interval of 0.1
* The normal mean is randomly chosen from a range between -0.5 and 0.9 with step interval of 0.1
* The normal standard deviation is randomly chosen from a range between 0.3 and 0.9 with step interval of 0.1

The mixture model is performed for a different number of random components, from 3 to 5. For each number of components, the EM is initialized and performed multiple times; then the best model (lowest BIC) is chosen as representative of the random method with such number of components and plotted in the figure. The criterion for the best model is 

Hybrid initialization
---------------------

The EM algorithm is initialized with the same components previously guessed from the *K*:sub:`S` data with the addition of a random lognormal component, based on the ranges used in the random initialization. As for the random initialization, the EM is initialized and performed multiple times and the best model is selected as representative of the hybrid method and plotted in the figure. The criterion for the best model is the lowest BIC.


Model evaluation
----------------

After having run the EM with all the three methods, the model with lowest BIC is considered the best one and plotted in a separate figure. The others are compared to it by the difference in their BIC scores (delta BIC).

Among the output files (for a complete list see section :ref:`output_files`), the ELMM produces a tabular (TSV) file and a text file where to store parameters and fitting results: 

* ``elmm_species_parameters.tsv``:

    * The type of model initialization is stored in column 1 according to a numerical code (1: data-driven, 2: data-driven plus a random lognormal component, 3: random initialization with exponential component and one lognormal component, 4: random initialization with exponential component and two lognormal components; higher numbers feature increasing number of lognormal components).
    * The initialization round is stored in column 2. By default each model type (execpt type 1) is initialized and fitted 10 times, so this column shows numbers from 1 to 10. 
    * The BIC and loglikelihood scores for the fitted model are stored in columns 3 and 4.
    * The number of EM iterations needed to reach convergence is stored in column 5. If greater than 300, the convergence is not reached and the cell shows *NA*.
    * The fitted exponential rate parameter and its component weight are stored in columns 6 and 7.
    * The mean, standard deviation and weight of the fitted Normal components used to define the correspondent lognormal components are stored in columns 8 to 10. When there are multiple lognormal components, the data for each of them are stored in a separate rows (the number of rows is thus equal to the number of lognormal components).

    .. figure:: _images/elmm.png
        :align: center
        :width: 800

        This file section shows the result for the first initalization of model 5: each row stores the same data for the exponential component plus the data for one of the three lognromal components.

* ``elmm_species_parameters.txt`` reports the results in a more descriptive and easy-to-read layout.


.. _`lmm`:

Lognormal mixture model
+++++++++++++++++++++++

The lognormal-only mixture modeling (LMM) uses only lognormal components and works by fitting Gaussians on the log-transformed *K*:sub:`S` distribution. The absence of the exponential component to model SSDs makes it less appropriate for paranome distributions, while this doesn't affect its application on anchor *K*:sub:`S` distributions. By default this method is turned off and can be switched on in the expert configuration file through ``extra_paralogs_analyses_methods``.


Among the output files (for a complete list see section :ref:`output_files`), the LMM produces tabular (TSV) files and text files where to store parameters and fitting results:

* ``lmm_species_parameters_paranome.tsv`` and ``lmm_species_parameters_anchors.tsv``:

    * The Model type is stored in column 1 according to a numerical code (1: one lognormal component, 2: two lognormal components, 3: three lognormal components; and so on).
    * The BIC and loglikelihood scores for the fitted model are stored in columns 2 and 3.
    * The number of EM iterations needed to reach convergence is stored in column 4. If greater than 300, the convergence is not reached and the cell shows *NA*.
    * The mean, standard deviation and weight of the fitted Normal components used to define the correspondent lognormal components are stored in columns 5 to 7. When there are multiple lognormal components, the data for each of them are stored in a separate rows (the number of rows is thus equal to the number of lognormal components).

    .. figure:: _images/lmm.png
        :align: center
        :width: 800

        This file section shows the result for model 5: each row stores the data for one of the five lognromal components.

* ``lmm_species_parameters_paranome.txt`` and ``lmm_species_parameters_anchors.txt`` collect the model results in a more descriptive and easy-to-read layout.
