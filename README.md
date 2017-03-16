# QGIS Scipy Clustering

This plugin implements point clustering in scipy and add a label integer field to the feature class for the clustered data. Both hierarchical and k-means clustering are implemented.

This is a Processing plugin (activated automatically) and can be found in the processing toolbox.

Please note that there are memory limitations in hierarchical clustering. The space required to create the clusters is O(n^2), which means that larger datasets will run out of memory fast. A plugin setting in the Processing options that sets the upper limit of points to process, by default set at 10,000. K-means is much more forgiving regarding memory, so the limit is not enforced in those algorithms.

All credit to the scipy team for the original implementation of the cluster
algorithms.

Jones E, Oliphant E, Peterson P, et al. SciPy: Open Source Scientific
Tools for Python, 2001-, http://www.scipy.org/ [Online].

seagull by Lane F. Kinkade from the Noun Project
https://thenounproject.com/term/seagull/166081
