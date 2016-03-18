#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Docstring
"""

__author__ = "henry"

# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load DifferentialPrivacy class from file DifferentialPrivacy.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .scipy_point_clustering import ScipyPointClusteringPlugin
    return ScipyPointClusteringPlugin()
