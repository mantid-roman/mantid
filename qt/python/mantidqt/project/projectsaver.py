# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
from __future__ import (absolute_import, division, print_function, unicode_literals)

from json import dump
import os

from mantidqt.project.workspacesaver import WorkspaceSaver
from mantidqt.project.plotssaver import PlotsSaver
from mantid import logger


class ProjectSaver(object):
    def __init__(self, project_file_ext):
        self.project_file_ext = project_file_ext

    def save_project(self, directory, workspace_to_save=None, interfaces_to_save=None, plots_to_save=None):
        """
        :param directory: String;
        :param workspace_to_save: List; of Strings that will have workspace names in it
        :param interfaces_to_save: List; of Strings that will have interface tags in it
        :param plots_to_save: List; of String that will have plot identifiers in it
        :return:
        """
        # Check if the directory doesn't exist
        if directory is None:
            logger.warning("Can not save to empty directory")
            return

        # Check this isn't saving a blank project file
        if workspace_to_save is None and interfaces_to_save is None and plots_to_save is None:
            logger.warning("Can not save an empty project")
            return

        # Save workspaces to that location
        workspace_saver = WorkspaceSaver(directory=directory)
        workspace_saver.save_workspaces(workspaces_to_save=workspace_to_save)

        # Generate plots
        plots_to_save_list = PlotsSaver().save_plots(plots_to_save)

        # Get interface details in dicts
        interfaces_to_save_dict = {}

        # Pass dicts to Project Writer
        writer = ProjectWriter(interfaces_to_save=interfaces_to_save_dict,
                               workspace_names=workspace_saver.get_output_list(),
                               plots_to_save=plots_to_save_list,
                               save_location=directory,
                               project_file_ext=self.project_file_ext)
        writer.write_out()


class ProjectWriter(object):
    def __init__(self, interfaces_to_save, workspace_names, plots_to_save, save_location, project_file_ext):
        self.interfaces_to_save = interfaces_to_save
        self.workspace_names = workspace_names
        self.directory = save_location
        self.project_file_ext = project_file_ext
        self.plots_to_save = plots_to_save

    def write_out(self):
        # Get the JSON string versions
        workspace_interface_dict = {"workspaces": self.workspace_names, "interfaces": self.interfaces_to_save,
                                    "plots": self.plots_to_save}

        # Open file and save the string to it alongside the workspace_names
        file_name = self.directory + '/' + (os.path.basename(self.directory) + self.project_file_ext)
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)
        f = open(file_name, 'w+')
        dump(obj=workspace_interface_dict, fp=f, indent=2)