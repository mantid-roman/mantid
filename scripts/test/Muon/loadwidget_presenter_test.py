import sys
import os
import time

from Muon.GUI.MuonAnalysis.loadwidget.load_widget_model import LoadWidgetModel
from Muon.GUI.MuonAnalysis.loadwidget.load_widget_view import LoadWidgetView
from Muon.GUI.MuonAnalysis.loadwidget.load_widget_presenter import LoadWidgetPresenter

from Muon.GUI.MuonAnalysis.loadrun.load_run_model_multithreading import LoadRunWidgetModel
from Muon.GUI.MuonAnalysis.loadrun.load_run_presenter import LoadRunWidgetPresenter
from Muon.GUI.MuonAnalysis.loadrun.load_run_view import LoadRunWidgetView

from Muon.GUI.MuonAnalysis.loadfile.load_file_model_multithreading import BrowseFileWidgetModel
from Muon.GUI.MuonAnalysis.loadfile.load_file_presenter import BrowseFileWidgetPresenter
from Muon.GUI.MuonAnalysis.loadfile.load_file_view import BrowseFileWidgetView

from Muon.GUI.Common.muon_load_data import MuonLoadData
import Muon.GUI.Common.muon_file_utils as fileUtils

# TODO : Fix these tests so that they finish with exit code 0 (rather than hanging)

import unittest

if sys.version_info.major == 3:
    from unittest import mock
else:
    import mock

from qtpy import QtWidgets


class LoadRunWidgetPresenterTest(unittest.TestCase):
    class Runner:
        QT_APP = QtWidgets.QApplication([])

        def __init__(self, thread):
            if thread:
                self._thread = thread
                self._thread.finished.connect(self.finished)
                if self._thread.is_running():
                    self.QT_APP.exec_()

        def finished(self):
            self.QT_APP.processEvents()
            self.QT_APP.exit(0)

    def setUp(self):
        self.data = MuonLoadData()
        self.load_file_view = BrowseFileWidgetView()
        self.load_run_view = LoadRunWidgetView()
        self.load_file_model = BrowseFileWidgetModel(self.data)
        self.load_run_model = LoadRunWidgetModel(self.data)

        self.presenter = LoadWidgetPresenter(
            LoadWidgetView(load_file_view=self.load_file_view, load_run_view=self.load_run_view),
            LoadWidgetModel(self.data))
        self.presenter.set_load_file_widget(BrowseFileWidgetPresenter(self.load_file_view, self.load_file_model))
        self.presenter.set_load_run_widget(LoadRunWidgetPresenter(self.load_run_view, self.load_run_model))

        self.mock_loading_from_browse([1], "C:\dir1\dir2\dir3\EMU0001234.nxs", 1234)
        fileUtils.get_current_run_filename = mock.Mock(return_value="EMU0001234.nxs")

    def mock_loading_from_browse(self, workspace, filename, run):
        self.load_file_view.show_file_browser_and_return_selection = mock.Mock(
            return_value=[filename])
        self.load_file_model.load_workspace_from_filename = mock.Mock(
            return_value=(filename, workspace, run))

    def mock_loading_from_current_run(self, workspace, filename, run):
        self.load_run_model.load_workspace_from_filename = mock.Mock(
            return_value=(workspace, filename, run))

    def mock_user_input_single_run(self, workspace, filename, run):
        self.load_run_view.get_run_edit_text = mock.Mock(return_value=str(run))
        self.load_run_model.load_workspace_from_filename = mock.Mock(
            return_value=(workspace, filename, run))

    def mock_user_input_single_file(self, workspace, filename, run):
        self.load_file_model.load_workspace_from_filename = mock.Mock(
            return_value=(filename, workspace, run))
        self.load_file_view.get_file_edit_text = mock.Mock(
            return_value=filename)

    def mock_disabling_buttons_in_run_and_file_widget(self):
        self.load_run_view.disable_load_buttons = mock.Mock()
        self.load_run_view.enable_load_buttons = mock.Mock()
        self.load_file_view.disable_load_buttons = mock.Mock()
        self.load_file_view.enable_load_buttons = mock.Mock()

    # ------------------------------------------------------------------------------------------------------------------
    # TESTS : The interface should correctly load/display data via the run widget (load current run and user input run)
    # and via the file widget (browse button or user enters file)
    # ------------------------------------------------------------------------------------------------------------------

    def test_that_loading_single_file_via_browse_sets_model_data(self):
        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.presenter._model.filenames, ["C:\dir1\dir2\dir3\EMU0001234.nxs"])
        self.assertEqual(self.presenter._model.workspaces, [[1]])
        self.assertEqual(self.presenter._model.runs, [1234])

    def test_that_loading_single_file_via_browse_sets_run_and_file_view(self):
        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.get_file_edit_text(), "C:\dir1\dir2\dir3\EMU0001234.nxs")
        self.assertEqual(self.load_run_view.get_run_edit_text(), "1234")

    def test_that_loading_via_browse_disables_all_load_buttons(self):
        self.mock_disabling_buttons_in_run_and_file_widget()

        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.disable_load_buttons.call_count,
                         self.load_file_view.enable_load_buttons.call_count)
        self.assertEqual(self.load_run_view.disable_load_buttons.call_count,
                         self.load_run_view.enable_load_buttons.call_count)

    def test_that_loading_single_file_via_user_input_file_is_loaded_into_model_correctly(self):
        self.mock_user_input_single_file([1], "C:\dir1\dir2\dir3\EMU0001111.nxs", 1111)

        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.presenter._model.filenames, ["C:\dir1\dir2\dir3\EMU0001111.nxs"])
        self.assertEqual(self.presenter._model.workspaces, [[1]])
        self.assertEqual(self.presenter._model.runs, [1111])

    def test_that_loading_single_file_via_user_input_file_is_displayed_correctly_in_the_interface(self):
        self.mock_user_input_single_file([1], "C:\dir1\dir2\dir3\EMU0001111.nxs", 1111)

        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.get_file_edit_text(), "C:\dir1\dir2\dir3\EMU0001111.nxs")
        self.assertEqual(self.load_run_view.get_run_edit_text(), "1111")

    def test_that_loading_via_user_input_file_disables_all_load_buttons(self):
        self.mock_user_input_single_file([1], "C:\dir1\dir2\dir3\EMU0001111.nxs", 1111)

        self.mock_disabling_buttons_in_run_and_file_widget()

        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.disable_load_buttons.call_count,
                         self.load_file_view.enable_load_buttons.call_count)
        self.assertEqual(self.load_run_view.disable_load_buttons.call_count,
                         self.load_run_view.enable_load_buttons.call_count)

    def test_that_loading_via_load_current_run_is_loaded_into_model_correctly(self):
        self.mock_loading_from_current_run([1], "\\\\EMU\\data\\EMU0083420.nxs", 83420)

        self.presenter.load_run_widget.handle_load_current_run()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.presenter._model.filenames, ["\\\\EMU\\data\\EMU0083420.nxs"])
        self.assertEqual(self.presenter._model.workspaces, [[1]])
        self.assertEqual(self.presenter._model.runs, [83420])

    def test_that_loading_via_load_current_run_is_displayed_correctly_in_the_interface(self):
        self.mock_loading_from_current_run([1], "\\\\EMU\\data\\EMU0083420.nxs", 83420)

        self.presenter.load_run_widget.handle_load_current_run()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.get_file_edit_text(), "\\\\EMU\\data\\EMU0083420.nxs")
        self.assertEqual(self.load_run_view.get_run_edit_text(), "83420")

    def test_that_loading_via_load_current_run_disables_all_load_buttons(self):
        self.mock_loading_from_current_run([1], "\\\\EMU\\data\\EMU0083420.nxs", 83420)
        self.mock_disabling_buttons_in_run_and_file_widget()

        self.presenter.load_run_widget.handle_load_current_run()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.disable_load_buttons.call_count,
                         self.load_file_view.enable_load_buttons.call_count)
        self.assertEqual(self.load_run_view.disable_load_buttons.call_count,
                         self.load_run_view.enable_load_buttons.call_count)

    def test_that_user_input_run_is_loaded_into_model_correctly(self):
        self.mock_user_input_single_run([1], "EMU00012345.nxs", 12345)

        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.presenter._model.filenames, ["EMU00012345.nxs"])
        self.assertEqual(self.presenter._model.workspaces, [[1]])
        self.assertEqual(self.presenter._model.runs, [12345])

    def test_that_user_input_run_is_displayed_correctly_in_the_interface(self):
        self.mock_user_input_single_run([1], "EMU00012345.nxs", 12345)

        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.get_file_edit_text(), "EMU00012345.nxs")
        self.assertEqual(self.load_run_view.get_run_edit_text(), "12345")

    def test_that_loading_via_user_input_run_disables_all_load_buttons(self):
        self.mock_user_input_single_run([1], "EMU00012345.nxs", 12345)
        self.mock_disabling_buttons_in_run_and_file_widget()

        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.disable_load_buttons.call_count,
                         self.load_file_view.enable_load_buttons.call_count)
        self.assertEqual(self.load_run_view.disable_load_buttons.call_count,
                         self.load_run_view.enable_load_buttons.call_count)


class LoadRunWidgetPresenterLoadFailTest(unittest.TestCase):
    class Runner:
        QT_APP = QtWidgets.QApplication([])

        def __init__(self, thread):
            if thread:
                self._thread = thread
                self._thread.finished.connect(self.finished)
                if self._thread.is_running():
                    self.QT_APP.exec_()

        def finished(self):
            self.QT_APP.processEvents()
            self.QT_APP.exit(0)

    def setUp(self):
        self.data = MuonLoadData()
        self.load_file_view = BrowseFileWidgetView()
        self.load_run_view = LoadRunWidgetView()
        self.load_file_model = BrowseFileWidgetModel(self.data)
        self.load_run_model = LoadRunWidgetModel(self.data)

        self.presenter = LoadWidgetPresenter(
            LoadWidgetView(load_file_view=self.load_file_view, load_run_view=self.load_run_view),
            LoadWidgetModel(self.data))
        self.presenter.set_load_file_widget(BrowseFileWidgetPresenter(self.load_file_view, self.load_file_model))
        self.presenter.set_load_run_widget(LoadRunWidgetPresenter(self.load_run_view, self.load_run_model))

        self.load_file_view.show_file_browser_and_return_selection = mock.Mock(
            return_value=["C:\\dir1\\EMU0001234.nxs"])
        self.load_file_model.load_workspace_from_filename = mock.Mock(return_value=("C:\\dir1\\EMU0001234.nxs",[1], 1234))

        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.mock_loading_to_throw()
        fileUtils.get_current_run_filename = mock.Mock(return_value="EMU0001234.nxs")

        self.load_file_view.warning_popup = mock.Mock()
        self.load_run_view.warning_popup = mock.Mock()

    def load_failure(self):
        raise ValueError("Error text")

    def mock_loading_to_throw(self):
        self.load_file_model.load_workspace_from_filename = mock.Mock(side_effect=self.load_failure)
        self.load_run_model.load_workspace_from_filename = mock.Mock(side_effect=self.load_failure)

    def mock_disabling_buttons_in_run_and_file_widget(self):
        self.load_run_view.disable_load_buttons = mock.Mock()
        self.load_run_view.enable_load_buttons = mock.Mock()
        self.load_file_view.disable_load_buttons = mock.Mock()
        self.load_file_view.enable_load_buttons = mock.Mock()

    def assert_model_unchanged(self):
        self.assertEqual(self.presenter._model.filenames, ["C:\\dir1\\EMU0001234.nxs"])
        self.assertEqual(self.presenter._model.workspaces, [[1]])
        self.assertEqual(self.presenter._model.runs, [1234])

    def assert_interface_unchanged(self):
        self.assertEqual(self.load_file_view.get_file_edit_text(), "C:\\dir1\\EMU0001234.nxs")
        self.assertEqual(self.load_run_view.get_run_edit_text(), "1234")

    # ------------------------------------------------------------------------------------------------------------------
    # TESTS : The interface should always revert to its previous state if a load fails from anywhere in the widget
    # and a warning should be shown to the user.
    # ------------------------------------------------------------------------------------------------------------------

    def test_that_if_load_fails_from_browse_that_model_and_interface_are_unchanged_from_previous_state(self):
        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assert_model_unchanged()
        self.assert_interface_unchanged()

    def test_that_if_load_fails_from_user_file_entry_that_model_and_interface_are_unchanged_from_previous_state(self):
        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assert_model_unchanged()
        self.assert_interface_unchanged()

    def test_that_if_load_fails_from_current_run_that_model_and_interface_are_unchanged_from_previous_state(self):
        self.presenter.load_run_widget.handle_load_current_run()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assert_model_unchanged()
        self.assert_interface_unchanged()

    def test_that_if_load_fails_from_user_run_entry_that_model_and_interface_are_unchanged_from_previous_state(self):
        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assert_model_unchanged()
        self.assert_interface_unchanged()

    @mock.patch("Muon.GUI.Common.thread_model.warning")
    def test_that_if_load_fails_from_browse_that_warning_is_displayed(self, mock_warning):
        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.warning_popup.call_count, 1)

    @mock.patch("Muon.GUI.Common.thread_model.warning")
    def test_that_if_load_fails_from_user_file_entry_that_warning_is_displayed(self, mock_warning):
        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assertEqual(self.load_file_view.warning_popup.call_count, 1)

    @mock.patch("Muon.GUI.Common.thread_model.warning")
    def test_that_if_load_fails_from_current_run_that_warning_is_displayed(self, mock_warning):
        self.presenter.load_run_widget.handle_load_current_run()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_run_view.warning_popup.call_count, 1)

    @mock.patch("Muon.GUI.Common.thread_model.warning")
    def test_that_if_load_fails_from_user_run_entry_that_warning_is_displayed(self, mock_warning):
        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assertEqual(self.load_run_view.warning_popup.call_count, 1)


class LoadRunWidgetPresenterMultipleFileTest(unittest.TestCase):
    class Runner:
        QT_APP = QtWidgets.QApplication([])

        def __init__(self, thread):
            if thread:
                self._thread = thread
                self._thread.finished.connect(self.finished)
                if self._thread.is_running():
                    self.QT_APP.exec_()

        def finished(self):
            self.QT_APP.processEvents()
            self.QT_APP.exit(0)

    def setUp(self):
        self.data = MuonLoadData()
        self.load_file_view = BrowseFileWidgetView()
        self.load_run_view = LoadRunWidgetView()
        self.load_file_model = BrowseFileWidgetModel(self.data)
        self.load_run_model = LoadRunWidgetModel(self.data)

        self.view = LoadWidgetView(load_file_view=self.load_file_view, load_run_view=self.load_run_view)
        self.presenter = LoadWidgetPresenter(
            self.view,
            LoadWidgetModel(self.data))
        self.presenter.set_load_file_widget(BrowseFileWidgetPresenter(self.load_file_view, self.load_file_model))
        self.presenter.set_load_run_widget(LoadRunWidgetPresenter(self.load_run_view, self.load_run_model))

        self.view.multiple_loading_check.setCheckState(1)
        self.presenter.handle_multiple_files_option_changed()

        self.runs = [1234 + i for i in range(10)]
        self.workspaces = [[run] for run in self.runs]
        self.filenames = ["C:\dir1\dir2\EMU000" + str(run) + ".nxs" for run in self.runs]

    def mock_loading_from_browse(self, workspaces, filenames, runs):
        self.load_file_view.show_file_browser_and_return_selection = mock.Mock(
            return_value=filenames)
        self.load_file_model.load_workspace_from_filename = mock.Mock(
            side_effect=iter(zip(filenames, workspaces, runs)))

    # def mock_loading_from_current_run(self, workspace, filename, run):
    #     self.load_run_model.load_workspace_from_filename = mock.Mock(
    #         return_value=(workspace, filename, run))
    #
    def mock_user_input_multiple_runs(self, user_input, workspaces, filenames, runs):
        self.load_run_view.get_run_edit_text = mock.Mock(return_value=user_input)
        self.load_run_model.load_workspace_from_filename = mock.Mock(
            side_effect=iter(zip(workspaces, filenames, runs)))

    def mock_user_input_multiple_files(self, workspaces, filenames, runs):
        self.load_file_view.get_file_edit_text = mock.Mock(
            return_value= ";".join(filenames))
        self.load_file_model.load_workspace_from_filename = mock.Mock(
            side_effect=iter(zip(filenames, workspaces, runs)))

    def assert_model_contains_correct_loaded_data(self):
        # use sorted due to threads finishing at different times
        self.assertEqual(sorted(self.presenter._model.filenames), sorted(self.filenames))
        self.assertEqual(sorted(self.presenter._model.workspaces), sorted(self.workspaces))
        self.assertEqual(sorted(self.presenter._model.runs), sorted(self.runs))

    def assert_interface_contains_correct_runs_and_files(self):
        self.assertEqual(self.load_file_view.get_file_edit_text(), ";".join(self.filenames))
        self.assertEqual(self.load_run_view.get_run_edit_text(), "1234-1243")


    def mock_disabling_buttons_in_run_and_file_widget(self):
        self.load_run_view.disable_load_buttons = mock.Mock()
        self.load_run_view.enable_load_buttons = mock.Mock()
        self.load_file_view.disable_load_buttons = mock.Mock()
        self.load_file_view.enable_load_buttons = mock.Mock()

    # ------------------------------------------------------------------------------------------------------------------
    # TESTS : The loading of multiple files is supported by the widget
    # ------------------------------------------------------------------------------------------------------------------

    def test_that_loading_multiple_files_via_browse_sets_model_and_interface_correctly(self):
        self.mock_loading_from_browse(self.workspaces, self.filenames, self.runs)

        self.presenter.load_file_widget.on_browse_button_clicked()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assert_model_contains_correct_loaded_data()
        self.assert_interface_contains_correct_runs_and_files()

    def test_that_loading_multiple_files_via_user_input_files_sets_model_and_interface_correctly(self):
        self.mock_user_input_multiple_files(self.workspaces, self.filenames, self.runs)

        self.presenter.load_file_widget.handle_file_changed_by_user()
        self.Runner(self.presenter.load_file_widget._model.thread_manager)

        self.assert_model_contains_correct_loaded_data()
        self.assert_interface_contains_correct_runs_and_files()

    def test_that_loading_multiple_files_via_user_input_runs_sets_model_and_interface_correctly(self):
        self.mock_user_input_multiple_runs("1234-1243", self.workspaces, self.filenames, self.runs)

        self.presenter.load_run_widget.handle_run_changed_by_user()
        self.Runner(self.presenter.load_run_widget._model.thread_manager)

        self.assert_model_contains_correct_loaded_data()
        self.assert_interface_contains_correct_runs_and_files()

    def test_that_loading_multiple_files_then_allows_loading_via_increment_and_decrement(self):
        pass


if __name__ == '__main__':
    unittest.main(buffer=False, verbosity=2)