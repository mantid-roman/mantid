# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
#
from __future__ import (absolute_import, unicode_literals)

# std imports
import os.path as osp

# 3rd party imports
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (QTabWidget, QToolButton, QVBoxLayout, QWidget)

# local imports
from mantidqt.widgets.codeeditor.interpreter import PythonFileInterpreter
from mantidqt.widgets.codeeditor.scriptcompatibility import (mantid_api_import_needed,
                                                             add_mantid_api_import)


NEW_TAB_TITLE = 'New'
MODIFIED_MARKER = '*'


def _tab_title_and_toolip(filename):
    """Create labels for the tab title and tooltip from a filename"""
    if filename is None:
        return NEW_TAB_TITLE, NEW_TAB_TITLE
    else:
        return osp.basename(filename), filename


class MultiPythonFileInterpreter(QWidget):
    """Provides a tabbed widget for editing multiple files"""

    def __init__(self, default_content=None, parent=None):
        super(MultiPythonFileInterpreter, self).__init__(parent)

        # attributes
        self.default_content = default_content
        self.prev_session_tabs = None
        self.whitespace_visible = False

        # widget setup
        self._tabs = self.create_tabwidget()
        layout = QVBoxLayout()
        layout.addWidget(self._tabs)
        self.setLayout(layout)
        layout.setContentsMargins(0, 0, 0, 0)

        # add a single editor by default
        self.append_new_editor()

    @property
    def editor_count(self):
        return self._tabs.count()

    @property
    def tab_filepaths(self):
        file_paths = []
        for idx in range(self.editor_count):
            file_path = self._tabs.widget(idx).filename
            if file_path:
                file_paths.append(file_path)
        return file_paths

    def append_new_editor(self, content=None, filename=None):
        if content is None:
            content = self.default_content
        interpreter = PythonFileInterpreter(content, filename=filename,
                                            parent=self._tabs)
        if self.whitespace_visible:
            interpreter.set_whitespace_visible()

        # monitor future modifications
        interpreter.sig_editor_modified.connect(self.mark_current_tab_modified)
        interpreter.sig_filename_modified.connect(self.on_filename_modified)

        tab_title, tab_tooltip = _tab_title_and_toolip(filename)
        tab_idx = self._tabs.addTab(interpreter, tab_title)
        self._tabs.setTabToolTip(tab_idx, tab_tooltip)
        self._tabs.setCurrentIndex(tab_idx)
        return tab_idx

    def abort_current(self):
        """Request that that the current execution be cancelled"""
        self.current_editor().abort()

    def close_all(self):
        """
        Close all tabs
        :return: True if all tabs are closed, False if cancelled
        """
        for idx in reversed(range(self.editor_count)):
            if not self.close_tab(idx):
                return False

        return True

    def close_tab(self, idx):
        """
        Close the tab at the given index.
        :param idx: The tab index
        :return: True if tab is to be closed, False if cancelled
        """
        if idx >= self.editor_count:
            return True
        # Make the current tab active so that it is clear what you
        # are being prompted to save
        self._tabs.setCurrentIndex(idx)
        if self.current_editor().confirm_close():
            widget = self._tabs.widget(idx)
            # note: this does not close the widget, that is why we manually close it
            self._tabs.removeTab(idx)
            widget.close()
        else:
            return False

        # we never want an empty widget
        if self.editor_count == 0:
            self.append_new_editor()

        return True

    def create_tabwidget(self):
        """Create a new QTabWidget with a button to add new tabs"""
        tabs = QTabWidget(self)
        tabs.setMovable(True)
        tabs.setTabsClosable(True)
        # create a button to add new tabs
        plus_btn = QToolButton(tabs)
        plus_btn.setText('+')
        plus_btn.clicked.connect(self.plus_button_clicked)
        tabs.setCornerWidget(plus_btn, Qt.TopLeftCorner)
        tabs.tabCloseRequested.connect(self.close_tab)
        return tabs

    def current_editor(self):
        return self._tabs.currentWidget()

    def editor_at(self, idx):
        """Return the editor at the given index. Must be in range"""
        return self._tabs.widget(idx)

    def execute_current(self):
        """Execute content of the current file. If a selection is active
        then only this portion of code is executed"""
        self.current_editor().execute_async()

    def mark_current_tab_modified(self, modified):
        """Update the current tab title to indicate that the
        content has been modified"""
        self.mark_tab_modified(self._tabs.currentIndex(), modified)

    def mark_tab_modified(self, idx, modified):
        """Update the tab title to indicate that the
        content has been modified or not"""
        title_cur = self._tabs.tabText(idx)
        if modified:
            if not title_cur.endswith(MODIFIED_MARKER):
                title_new = title_cur + MODIFIED_MARKER
            else:
                title_new = title_cur
        else:
            if title_cur.endswith(MODIFIED_MARKER):
                title_new = title_cur.rstrip('*')
            else:
                title_new = title_cur
        self._tabs.setTabText(idx, title_new)

    def on_filename_modified(self, filename):
        title, tooltip = _tab_title_and_toolip(filename)
        idx_cur = self._tabs.currentIndex()
        self._tabs.setTabText(idx_cur, title)
        self._tabs.setTabToolTip(idx_cur, tooltip)

    def open_file_in_new_tab(self, filepath, startup=False):
        """Open the existing file in a new tab in the editor

        :param filepath: A path to an existing file
        :param startup: Flag for if function is being called on startup
        """
        with open(filepath, 'r') as code_file:
            content = code_file.read()

        self.append_new_editor(content=content, filename=filepath)
        if startup is False and mantid_api_import_needed(content) is True:
            add_mantid_api_import(self.current_editor().editor, content)

    def open_files_in_new_tabs(self, filepaths):
        for filepath in filepaths:
            self.open_file_in_new_tab(filepath)

    def plus_button_clicked(self, _):
        """Add a new tab when the plus button is clicked"""
        self.append_new_editor()

    def restore_session_tabs(self):
        if self.prev_session_tabs is not None:
            try:
                self.open_files_in_new_tabs(self.prev_session_tabs)
            except IOError:
                pass
            self.close_tab(0)  # close default empty script

    def save_current_file(self):
        """Save the current file"""
        self.current_editor().save()

    def spaces_to_tabs_current(self):
        self.current_editor().replace_spaces_with_tabs()

    def tabs_to_spaces_current(self):
        self.current_editor().replace_tabs_with_spaces()

    def toggle_comment_current(self):
        self.current_editor().toggle_comment()

    def toggle_find_replace_dialog(self):
        self.current_editor().show_find_replace_dialog()

    def toggle_whitespace_visible_all(self):
        if self.whitespace_visible:
            for idx in range(self.editor_count):
                self.editor_at(idx).set_whitespace_invisible()
            self.whitespace_visible = False
        else:
            for idx in range(self.editor_count):
                self.editor_at(idx).set_whitespace_visible()
            self.whitespace_visible = True
