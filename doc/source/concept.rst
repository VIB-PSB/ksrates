*K*:sub:`S` rate-adjustment
***************************

Substitution rate-adjustment strategy in a nutshell
===================================================

``ksrates`` is a package for substitution rate-adjustment in mixed ortholog and paralog *K*:sub:`S` distributions.

Mixed *K*:sub:`S` distributions are one of the approaches applied to detect whole-genome duplications (WGDs) and to locate them in a phylogeny. A mixed plot is composed of ortholog *K*:sub:`S` distributions - representing divergence events - overlapped onto paralog *K*:sub:`S` distributions - representing the duplication history of a species genome. The relative positions of the ortholog peaks and the WGDs peaks are informative about the order of the depicted evolutionary events, allowing to place the occurrence of a WGDs in a specific branch of the evolutionary history of the species.

The reliability of a mixed plot can be jeopardized in case of (remarkable) substitution rate differences between the involved species. In fact, since the *K*:sub:`S` value of a homolog pair depends on the substitution rate of the species, different distributions end up to be built on different *K*:sub:`S` scales. A direct overlap of distributions is therefore likely to lead to unreliable interpretations.

The *K*:sub:`S` rate-adjustment package offers an adjustment procedure that brings all the distributions to a common *K*:sub:`S` scale by compensating for the substitution rate differences relatively to one "main" species. 
The rate-adjusted mixed plot obtained through ``ksrates`` is composed of a) a single paralog distribution coming from the main species and b) one or more ortholog distributions between the main species and the another species. The analysis is thus focused on the genome duplication history of the main species in the context of its evolutionary history with the other species. 

The rate-adjustment is applied to all the ortholog distributions. For each ortholog distribution, principles from the relative rate test (RRT) are used to detect the relative rates between the main species and the other species. During the rate-adjustment, the ortholog *K*:sub:`S`  peak is re-encoded as twice the relative rate of the main species, so that the age of the ortholog distribution is adapted to the *K*:sub:`S` scale of the paralog distribution. At the end, all ortholog distributions are seen from the perspective of the main species rate.
The rate-adjustment generates horizontal shifts of the ortholog distribution peak towards left if the main species is slower than the other species, or towards right if it is faster. The new disposition of the divergence events can lead to a different and more reliable interpretation of WGD placement or of the order of the divergences themselves.
For more details about the rate-adjustment strategy, see [...].


.. _`explained_example`:

Explained example
=================

This example studies the phylogenetic placement of WGD signatures present in oil palm (*Elaeis guineensis*) paralog distribution. The rate-adjustment pipeline needs a input phylogenetic tree and the sequence data of all involved species. The minimum input tree is composed by the focal species (palm), another species (rice) and their outgroup (asparagus): ``((palm, rice)), asparagus)``. 
..  The mixed plot will show the palm paralog distribution overlapped with the rate-adjusted ortholog distributions involving palm and the other species in the input tree.

From the perspective of palm history there are two divergence events (i.e. ortholog distributions) in this tree, namely palm-rice and palm-asparagus. The pipeline breaks down the tree into *trios* composed by the species pair of a ortholog distribution and an outgroup used for its rate-adjustment. The example tree gives only one trio, "palm, rice, asparagus", where palm-rice divergence is rate-adjusted with outgroup asparagus. Palm-asparagus divergence has instead no outgroup in this tree and will be ignored; to avoid this, add another outgroup to the phylogeny, e.g. ``(((palm, rice), asparagus), spirodela)``. The user can also decide to perform multiple rate-adjustments for a divergence if the tree structure allows it: for example in this latter tree palm-rice can be rate-adjusted both with asparagus and spirodela (*Spirodela polyrhiza*).

Further on, the pipeline breaks down the trios into the three possible species pairs they are composed of, which in this case are palm-rice, palm-asparagus and rice-asparagus. ``wgd`` package then estimates the ortholog *K*:sub:`S` distribution for each of them. The ortholog distributions are simplified to a vertical line centered on their peak value (Figure 1).

.. figure:: _images/ortholog_distribution_peak.svg
    :align: center
    :width: 350

    The ortholog distribution for palm and rice is approximated to its mode (1.53 *K*:sub:`S`).
    
The RRT uses the *K*:sub:`S` values of the three ortholog peaks to compute the relative rates of the divergent pair: palm has a relative rate of about 0.36 while rice of 1.17, therefore palm accumulates substitution much more slowly than rice. Lastly, the rate-adjustment reinterprets the ortholog *K*:sub:`S` peak of palm-rice by encoding it as twice the relative rate of palm (*K*:sub:`S`' = 0.73). The ortholog peak has therefore been largely shifted to the left from 1.53 to 0.73 *K*:sub:`S` (Figure 2), and it is now adapted to the slow scale of palm paralog distribution. The shift has important consequences in the interpretation of the mixed plot concerning the older WGD signal around 0.9 *K*:sub:`S`.

.. figure:: _images/mixed_palm_corrected.svg
    :align: center
    :width: 800

    The ortholog distribution peak (red line) has been shifted towards left after rate-adjustment, as highlighted by the red arrows starting from the original position and pointing at the new rate-adjusted position. 

