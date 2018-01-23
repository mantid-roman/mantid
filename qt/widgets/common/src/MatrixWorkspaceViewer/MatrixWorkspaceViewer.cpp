#include "MantidQtWidgets/Common/MatrixWorkspaceViewer/MatrixWorkspaceViewer.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/AnalysisDataService.h"
#include <MantidQtWidgets/Common/pixmaps.h>

#include <QMessageBox>
#include <QHeaderView>
#include <QVBoxLayout>
#include <QScrollBar>
#include <QTableView>

#include <iostream>

using namespace Mantid::API;
using namespace MantidQt::MantidWidgets;

MatrixWorkspaceViewer::MatrixWorkspaceViewer(QString wsName, QWidget *parent) : 
  QFrame(parent), WorkspaceObserver(), m_layout(nullptr)
{
  auto ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(wsName.toStdString());
  setWorkspace(ws);
}

void MatrixWorkspaceViewer::setup(
  Mantid::API::MatrixWorkspace_sptr ws,
  int start,
  int end)
{
  if (!ws) {
    QMessageBox::critical(0, "WorkspaceMatrixModel error",
                          "2D workspace expected.");
    m_rows = 0;
    m_cols = 0;
    m_startRow = 0;
    m_endRow = 0;
    return;
  }

  m_workspace = ws;
  m_workspaceTotalHist = static_cast<int>(ws->getNumberHistograms());
  m_startRow = (start < 0 || start >= m_workspaceTotalHist) ? 0 : start;
  m_endRow = (end < 0 || end >= m_workspaceTotalHist || end < start)
                 ? m_workspaceTotalHist - 1
                 : end;
  m_rows = m_endRow - m_startRow + 1;
  try {
    // let the workspace do its thing
    m_cols = static_cast<int>(ws->blocksize());
  } catch (std::length_error &) {
    // otherwise get the maximum
    m_cols = static_cast<int>(ws->y(0).size());
    for (int i = 0; i < m_workspaceTotalHist; ++i) {
      m_cols = std::max(m_cols, static_cast<int>(ws->y(i).size()));
    }
  }
  if (ws->isHistogramData())
    m_histogram = true;
  connect(this, SIGNAL(needsUpdating()), this, SLOT(repaintAll()));

  m_bk_color = QColor(128, 255, 255);
  m_matrix_icon = API::getQPixmap("mantid_matrix_xpm");
  m_column_width = 100;
}


void MatrixWorkspaceViewer::setWorkspace(Mantid::API::MatrixWorkspace_sptr ws, int start, int end) {
  m_workspace = ws;

  setup(ws, start, end);

  m_modelY = new MantidMatrixModel(this, ws.get(), m_rows, m_cols, m_startRow,
                                   MantidMatrixModel::Y);
  m_table_viewY = new QTableView();
  connectTableView(m_table_viewY, m_modelY);
  //setColumnsWidth(0, MantidPreferences::MantidMatrixColumnWidthY());
  //setNumberFormat(0, MantidPreferences::MantidMatrixNumberFormatY(),
  //                MantidPreferences::MantidMatrixNumberPrecisionY());

  m_modelX = new MantidMatrixModel(this, ws.get(), m_rows, m_cols, m_startRow,
                                   MantidMatrixModel::X);
  m_table_viewX = new QTableView();
  connectTableView(m_table_viewX, m_modelX);
  //setColumnsWidth(1, MantidPreferences::MantidMatrixColumnWidthX());
  //setNumberFormat(1, MantidPreferences::MantidMatrixNumberFormatX(),
  //                MantidPreferences::MantidMatrixNumberPrecisionX());

  m_modelE = new MantidMatrixModel(this, ws.get(), m_rows, m_cols, m_startRow,
                                   MantidMatrixModel::E);
  m_table_viewE = new QTableView();
  connectTableView(m_table_viewE, m_modelE);
  //setColumnsWidth(2, MantidPreferences::MantidMatrixColumnWidthE());
  //setNumberFormat(2, MantidPreferences::MantidMatrixNumberFormatE(),
  //                MantidPreferences::MantidMatrixNumberPrecisionE());

  m_YTabLabel = QString("Y values");
  m_XTabLabel = QString("X values");
  m_ETabLabel = QString("Errors");

  if (!m_layout) {
    m_layout = new QVBoxLayout();
    setLayout(m_layout);
  } else {
    m_layout->takeAt(0);
    m_tabs->close();
  }

  m_tabs = new QTabWidget(this);
  m_tabs->insertTab(0, m_table_viewY, m_YTabLabel);
  m_tabs->insertTab(1, m_table_viewX, m_XTabLabel);
  m_tabs->insertTab(2, m_table_viewE, m_ETabLabel);

  m_layout->addWidget(m_tabs);

  // for synchronizing the views
  // index is zero for the defualt view
  m_PrevIndex = 0;
  // install event filter on  these objects
  m_table_viewY->installEventFilter(this);
  m_table_viewX->installEventFilter(this);
  m_table_viewE->installEventFilter(this);

  connect(m_tabs, SIGNAL(currentChanged(int)), this, SLOT(viewChanged(int)));

  setGeometry(50, 50,
              qMin(5, numCols()) *
                      m_table_viewY->horizontalHeader()->sectionSize(0) +
                  55,
              (qMin(10, numRows()) + 1) *
                      m_table_viewY->verticalHeader()->sectionSize(0) +
                  100);

  // Add an extension for the DX component if required
  //if (ws->hasDx(0)) {
  //  addMantidMatrixTabExtension(MantidMatrixModel::DX);
  //}

  observeAfterReplace();
  observePreDelete();
  observeADSClear();

  connect(this, SIGNAL(needWorkspaceChange(Mantid::API::MatrixWorkspace_sptr)),
          this, SLOT(changeWorkspace(Mantid::API::MatrixWorkspace_sptr)));
  connect(this, SIGNAL(needToClose()), this, SLOT(closeMatrix()));

  connect(this, SIGNAL(closedWindow(MdiSubWindow *)), this,
          SLOT(selfClosed(MdiSubWindow *)));
}

void MatrixWorkspaceViewer::connectTableView(QTableView *view,
                                    MantidMatrixModel *model) {
  view->setSizePolicy(
      QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));
  view->setSelectionMode(QAbstractItemView::ExtendedSelection);
  view->setModel(model);
  view->setCornerButtonEnabled(false);
  view->setFocusPolicy(Qt::StrongFocus);

  QPalette pal = view->palette();
  pal.setColor(QPalette::Base, m_bk_color);
  view->setPalette(pal);

  // set header properties
  QHeaderView *hHeader = (QHeaderView *)view->horizontalHeader();
#if QT_VERSION < 0x050000
  hHeader->setMovable(false);
  hHeader->setResizeMode(QHeaderView::Interactive);
#else
  hHeader->setSectionsMovable(false);
  hHeader->setSectionResizeMode(QHeaderView::Interactive);
#endif
  hHeader->setDefaultSectionSize(m_column_width);

  view->resizeRowToContents(0);
  int row_height = view->rowHeight(0);

  QHeaderView *vHeader = (QHeaderView *)view->verticalHeader();
  vHeader->setDefaultSectionSize(row_height);
#if QT_VERSION < 0x050000
  vHeader->setResizeMode(QHeaderView::Fixed);
  vHeader->setMovable(false);
#else
  vHeader->setSectionResizeMode(QHeaderView::Fixed);
  vHeader->setSectionsMovable(false);
#endif
}

/** Called when switching between tabs
*  @param index :: The index of the new active tab
*/
void MatrixWorkspaceViewer::viewChanged(int index) {
  // get the previous view and selection model
  QTableView *prevView = (QTableView *)m_tabs->widget(m_PrevIndex);
  if (prevView) {
    QItemSelectionModel *oldSelModel = prevView->selectionModel();
    QItemSelectionModel *selModel = activeView()->selectionModel();
    // Copy the selection from the previous tab into the newly-activated one
    selModel->select(oldSelModel->selection(), QItemSelectionModel::Select);
    // Clear the selection on the now-hidden tab
    oldSelModel->clearSelection();

    m_PrevIndex = index;
    // get the previous tab scrollbar positions
    int hValue = prevView->horizontalScrollBar()->value();
    int vValue = prevView->verticalScrollBar()->value();
    // to synchronize the views
    // set  the previous view  scrollbar positions to current view
    activeView()->horizontalScrollBar()->setValue(hValue);
    activeView()->verticalScrollBar()->setValue(vValue);
  }
}

/**  Return the pointer to the active table view.
*/
QTableView *MatrixWorkspaceViewer::activeView() {
  auto currentIndex = m_tabs->currentIndex();
  switch (currentIndex) {
  case 0:
    return m_table_viewY;
  case 1:
    return m_table_viewX;
  case 2:
    return m_table_viewE;
  //default:
  //  return m_extensionRequest.getActiveView(intToModelType(currentIndex),
  //                                          m_extensions, m_table_viewY);
  }
}