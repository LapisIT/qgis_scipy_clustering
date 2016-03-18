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

import os
import sys
import inspect

from processing.core.Processing import Processing
from scipy_point_clustering_provider import ScipyPointClusteringProvider

cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]

if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class ScipyPointClusteringPlugin:

    def __init__(self):
        self.provider = ScipyPointClusteringProvider()

    def initGui(self):
        Processing.addProvider(self.provider)

    def unload(self):
        Processing.removeProvider(self.provider)
