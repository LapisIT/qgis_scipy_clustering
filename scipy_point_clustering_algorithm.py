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

from PyQt4.QtCore import QVariant
from qgis.core import QgsField, QgsFeature, QgsGeometry, QgsPoint, QgsFields

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.parameters import (
    ParameterVector, ParameterString, ParameterSelection, ParameterNumber,
    ParameterTableField
)
from processing.core.outputs import OutputVector
from processing.tools import dataobjects, vector
import numpy as np
import scipy.cluster.vq
import scipy.cluster.hierarchy
from scipy.spatial.distance import pdist, squareform

from scipy_point_clustering_utils import ScipyPointClusteringUtils


class HierarchicalClustering(GeoAlgorithm):
    """
    Implementation of hierarchical clustering from scipy.

    Works with point data and adds a label field to the original dataset to
    store the cluster.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT_LAYER = 'OUTPUT_LAYER'
    INPUT_LAYER = 'INPUT_LAYER'
    LINKAGE_METHOD = 'LINKAGE_METHOD'
    LINKAGE_METRIC = 'LINKAGE_METRIC'
    LABEL_FIELD = 'LABEL_FIELD'
    TOLERANCE = 'TOLERANCE'
    CRITERION = 'CRITERION'

    _linkage_methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']
    _linkage_metrics = ['euclidean', 'cityblock']
    _criterions = ['distance', 'inconsistent', 'maxclust', 'monocrit', 'maxclust_monocrit']

    def defineCharacteristics(self):
        """Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # The name that the user will see in the toolbox
        self.name = 'Hierarchical Clustering'

        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Vector'

        # We add the input vector layer. It can have any kind of geometry
        # It is a mandatory (not optional) one, hence the False argument
        self.addParameter(ParameterVector(self.INPUT_LAYER,
            self.tr('Input layer'), [ParameterVector.VECTOR_TYPE_POINT], False))

        self.addParameter(ParameterNumber(
            self.TOLERANCE, self.tr('Cluster tolerance'), minValue=0.0))

        self.addParameter(ParameterString(
            self.LABEL_FIELD, self.tr('Label field name'), 'label'
        ))

        self.addParameter(ParameterSelection(
            self.LINKAGE_METHOD, self.tr('Linkage method'),
            self._linkage_methods
        ))

        self.addParameter(ParameterSelection(
            self.LINKAGE_METRIC, self.tr('Linkage metric'),
            self._linkage_metrics
        ))

        self.addParameter(ParameterSelection(
            self.CRITERION, self.tr('Cluster criterion'),
            self._criterions
        ))

        # We add a vector layer as output
        self.addOutput(OutputVector(self.OUTPUT_LAYER,
            self.tr('Clustered features')))

    def processAlgorithm(self, progress):
        """Here is where the processing itself takes place."""

        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        inputFilename = self.getParameterValue(self.INPUT_LAYER)
        output = self.getOutputFromName(self.OUTPUT_LAYER)
        fieldName = self.getParameterValue(self.LABEL_FIELD)
        tolerance = float(self.getParameterValue(self.TOLERANCE))
        method = self._linkage_methods[self.getParameterValue(self.LINKAGE_METHOD)]
        metric = self._linkage_metrics[self.getParameterValue(self.LINKAGE_METRIC)]
        criterion = self._criterions[self.getParameterValue(self.CRITERION)]

        # Input layers vales are always a string with its location.
        # That string can be converted into a QGIS object (a
        # QgsVectorLayer in this case) using the
        # processing.getObjectFromUri() method.
        vectorLayer = dataobjects.getObjectFromUri(inputFilename)

        # And now we can process

        # Loop over the features to get the geometries and the associated
        # feature id
        features = vector.features(vectorLayer)
        points = []
        feature_ids = []
        for f in features:
            g = f.geometry()
            p = g.asPoint()
            points.append([p.x(), p.y()])
            feature_ids.append(f.id())

        # actually do the clustering
        points = np.array(points)

        y = scipy.cluster.hierarchy.fclusterdata(
            points,
            tolerance,
            criterion=criterion,
            method=method,
            metric=metric
        )

        labels = dict(zip(feature_ids, y))

        # Now write the features to the new dataset along with the label
        #label_idx = fields.fieldNameIndex(fieldName)

        provider = vectorLayer.dataProvider()
        fields = provider.fields()
        fields.append(QgsField(fieldName, QVariant.Int))
        writer = output.getVectorWriter(
            fields, provider.geometryType(), provider.crs())

        #assert False, labels

        out_feature = QgsFeature()
        out_feature.setFields(fields)

        features = vector.features(vectorLayer)
        for f in features:

            attributes = f.attributes()
            attributes.append(int(labels[f.id()]))

            geom = f.geometry()

            out_feature.setGeometry(geom)
            out_feature.setAttributes(attributes)

            writer.addFeature(out_feature)
        del writer

    def getIcon(self):
        """Get the icon.
        """
        return ScipyPointClusteringUtils.getIcon()


class KMeansClustering(GeoAlgorithm):
    """
    K-means clustering implementation
    """

    OUTPUT_LAYER = 'OUTPUT_LAYER'
    INPUT_LAYER = 'INPUT_LAYER'
    K = 'K'
    MINIT = 'MINIT'
    LABEL_FIELD = 'LABEL_FIELD'
    CENTROID_OUTPUT = 'CENTROID_OUTPUT'

    _minits = ['random', 'points']

    def defineCharacteristics(self):
        """Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        self.name = 'K-means clustering'

        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Vector'

        # We add the input vector layer. It can have any kind of geometry
        # It is a mandatory (not optional) one, hence the False argument
        self.addParameter(ParameterVector(
            self.INPUT_LAYER,
            self.tr('Input layer'),
            [ParameterVector.VECTOR_TYPE_POINT],
            False
        ))

        self.addParameter(ParameterNumber(
            self.K,
            self.tr('K (number of clusters'),
            minValue=2,
            default=3
        ))

        self.addParameter(ParameterString(
            self.LABEL_FIELD, self.tr('Label field name'), 'label'
        ))

        self.addParameter(ParameterSelection(
            self.MINIT,
            self.tr("Method for initialization"),
            self._minits
        ))

        # We add a vector layer as output
        self.addOutput(OutputVector(self.OUTPUT_LAYER,
                                    self.tr('Clustered features')))

        self.addOutput(OutputVector(self.CENTROID_OUTPUT,
                                    self.tr('Cluster centroids')))

    def processAlgorithm(self, progress):
        """Here is where the processing itself takes place."""

        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        inputFilename = self.getParameterValue(self.INPUT_LAYER)
        k = int(self.getParameterValue(self.K))
        minit = self._minits[self.getParameterValue(self.MINIT)]
        fieldName = self.getParameterValue(self.LABEL_FIELD)

        output = self.getOutputFromName(self.OUTPUT_LAYER)
        centroid_output = self.getOutputFromName(self.CENTROID_OUTPUT)

        # Input layers vales are always a string with its location.
        # That string can be converted into a QGIS object (a
        # QgsVectorLayer in this case) using the
        # processing.getObjectFromUri() method.
        vectorLayer = dataobjects.getObjectFromUri(inputFilename)

        # And now we can process

        # Loop over the features to get the geometries and the associated
        # feature id
        features = vector.features(vectorLayer)
        points = []
        feature_ids = []
        for f in features:
            g = f.geometry()
            p = g.asPoint()
            points.append([p.x(), p.y()])
            feature_ids.append(f.id())

        # actually do the clustering
        points = np.array(points)

        centroids, y = scipy.cluster.vq.kmeans2(points, k, minit=minit)

        labels = dict(zip(feature_ids, y))

        provider = vectorLayer.dataProvider()
        fields = provider.fields()
        fields.append(QgsField(fieldName, QVariant.Int))
        writer = output.getVectorWriter(
            fields, provider.geometryType(), provider.crs())

        # assert False, labels

        out_feature = QgsFeature()
        out_feature.setFields(fields)

        features = vector.features(vectorLayer)
        for f in features:
            attributes = f.attributes()
            attributes.append(int(labels[f.id()]))

            geom = f.geometry()

            out_feature.setGeometry(geom)
            out_feature.setAttributes(attributes)

            writer.addFeature(out_feature)
        del writer

        fields = QgsFields()
        fields.append(QgsField(fieldName, QVariant.Int))
        writer = centroid_output.getVectorWriter(
            fields, provider.geometryType(), provider.crs())

        out_feature = QgsFeature()
        out_feature.setFields(fields)

        for i, centroid in enumerate(centroids):

            geom = QgsPoint(*centroid)
            out_feature.setGeometry(QgsGeometry.fromPoint(geom))
            out_feature.setAttributes([i, ])
            writer.addFeature(out_feature)

        del writer


    def getIcon(self):
        """Get the icon.
        """
        return ScipyPointClusteringUtils.getIcon()


class HierarchicalClusteringByIdentifier(HierarchicalClustering):
    """
    Heirarchical clustering for features with some identifier.
    """

    IDENTIFIER_FIELD = 'IDENTIFIER_FIELD'

    def defineCharacteristics(self):
        """
        Call the parent and set attributes, then modify for this algorithm.
        """
        HierarchicalClustering.defineCharacteristics(self)

        self.name = 'Hierarchical Clustering by Identifier'

        self.parameters.insert(1, ParameterTableField(
            self.IDENTIFIER_FIELD, "Identifier field", self.INPUT_LAYER
        ))

    def processAlgorithm(self, progress):
        """Here is where the processing itself takes place."""

        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        inputFilename = self.getParameterValue(self.INPUT_LAYER)
        output = self.getOutputFromName(self.OUTPUT_LAYER)
        fieldName = self.getParameterValue(self.LABEL_FIELD)
        tolerance = float(self.getParameterValue(self.TOLERANCE))
        method = self._linkage_methods[
            self.getParameterValue(self.LINKAGE_METHOD)]
        metric = self._linkage_metrics[
            self.getParameterValue(self.LINKAGE_METRIC)]
        criterion = self._criterions[self.getParameterValue(self.CRITERION)]
        identifier_field = self.getParameterValue(self.IDENTIFIER_FIELD)

        # Input layers vales are always a string with its location.
        # That string can be converted into a QGIS object (a
        # QgsVectorLayer in this case) using the
        # processing.getObjectFromUri() method.
        vectorLayer = dataobjects.getObjectFromUri(inputFilename)

        # And now we can process

        # Loop over the features to get the geometries and the associated
        # feature id
        features = vector.features(vectorLayer)
        points = []
        feature_ids = []
        identifiers = []
        for f in features:
            g = f.geometry()
            p = g.asPoint()
            points.append([p.x(), p.y()])
            feature_ids.append(f.id())
            identifiers.append(f[identifier_field])

        # actually do the clustering
        points = np.array(points)
        identifiers = np.array(identifiers)

        # no we ensure that no matter how close the points are, the locatiosn of
        # the clusters is dependent on the label.
        distances = squareform(pdist(points, metric=metric))
        increase_offseets = np.array([
            identifiers != val for val in identifiers
        ])
        distances[increase_offseets] = np.inf
        distances = squareform(distances) # recompress the distance matrix

        links = scipy.cluster.hierarchy.linkage(
            distances, method=method, metric=metric)

        y = scipy.cluster.hierarchy.fcluster(
            links,
            tolerance,
            criterion=criterion
        )

        labels = dict(zip(feature_ids, y))

        # Now write the features to the new dataset along with the label
        # label_idx = fields.fieldNameIndex(fieldName)

        provider = vectorLayer.dataProvider()
        fields = provider.fields()
        fields.append(QgsField(fieldName, QVariant.Int))
        writer = output.getVectorWriter(
            fields, provider.geometryType(), provider.crs())

        # assert False, labels

        out_feature = QgsFeature()
        out_feature.setFields(fields)

        features = vector.features(vectorLayer)
        for f in features:
            attributes = f.attributes()
            attributes.append(int(labels[f.id()]))

            geom = f.geometry()

            out_feature.setGeometry(geom)
            out_feature.setAttributes(attributes)

            writer.addFeature(out_feature)
        del writer
