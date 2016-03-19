# -*- coding: utf-8 -*-

"""
/***************************************************************************
 ScipyPointClustering
                                 A QGIS plugin
 This plugin implements clustering for point data using the scipy module.
                              -------------------
        begin                : 2016-03-18
        copyright            : (C) 2016 by Henry Walshaw
        email                : henry.walshaw@spatialvision.com.au
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Henry Walshaw'
__date__ = '2016-03-18'
__copyright__ = '(C) 2016 by Henry Walshaw'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from processing.core.AlgorithmProvider import AlgorithmProvider
from processing.core.ProcessingConfig import Setting, ProcessingConfig
from scipy_point_clustering_algorithm import (
    HierarchicalClustering, KMeansClustering, HierarchicalClusteringByIdentifier
)
from scipy_point_clustering_utils import ScipyPointClusteringUtils


class ScipyPointClusteringProvider(AlgorithmProvider):


    def __init__(self):
        AlgorithmProvider.__init__(self)

        # Deactivate provider by default
        self.activate = True

        # Load algorithms
        self.alglist = [HierarchicalClustering(), KMeansClustering(), HierarchicalClusteringByIdentifier()]
        for alg in self.alglist:
            alg.provider = self

    def initializeSettings(self):
        """In this method we add settings needed to configure our
        provider.

        Do not forget to call the parent method, since it takes care
        or automatically adding a setting for activating or
        deactivating the algorithms in the provider.
        """
        AlgorithmProvider.initializeSettings(self)
        ProcessingConfig.addSetting(Setting(
            self.getDescription(),
            ScipyPointClusteringUtils.POINT_LIMIT,
            self.tr('Hierarchical clustering point limit'),
            10000
        ))

    def unload(self):
        """Setting should be removed here, so they do not appear anymore
        when the plugin is unloaded.
        """
        AlgorithmProvider.unload(self)
        ProcessingConfig.removeSetting(
            ScipyPointClusteringUtils.POINT_LIMIT
        )

    def getName(self):
        """This is the name that will appear on the toolbox group.

        It is also used to create the command line name of all the
        algorithms from this provider.
        """
        return 'Scipy Point Clustering'

    def getDescription(self):
        """This is the provired full name.
        """
        return 'Scipy Point Clustering'

    def getIcon(self):
        """Get the icon.
        """
        return ScipyPointClusteringUtils.getIcon()

    def _loadAlgorithms(self):
        """Here we fill the list of algorithms in self.algs.

        This method is called whenever the list of algorithms should
        be updated. If the list of algorithms can change (for instance,
        if it contains algorithms from user-defined scripts and a new
        script might have been added), you should create the list again
        here.

        In this case, since the list is always the same, we assign from
        the pre-made list. This assignment has to be done in this method
        even if the list does not change, since the self.algs list is
        cleared before calling this method.
        """
        self.algs = self.alglist
