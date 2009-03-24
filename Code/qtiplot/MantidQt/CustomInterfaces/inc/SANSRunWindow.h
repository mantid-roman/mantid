#ifndef MANTIDQTCUSTOMINTERFACES_SANSRUNWINDOW_H_
#define MANTIDQTCUSTOMINTERFACES_SANSRUNWINDOW_H_

//----------------------
// Includes
//----------------------
#include "MantidQtCustomInterfaces/ui_SANSRunWindow.h"
#include "MantidQtAPI/UserSubWindow.h"

#include <QMap>
#include <QHash>

//---------------------------
// Qt Forward Declarations
//---------------------------
class QLineEdit;

namespace MantidQt
{
namespace CustomInterfaces
{
//-----------------------------
// Forward declaration
//-----------------------------


class SANSRunWindow : public MantidQt::API::UserSubWindow
{
  Q_OBJECT

public:
  /// Default Constructor
  SANSRunWindow(QWidget *parent = 0);

  ~SANSRunWindow();

private:
  /// Initialize the layout
  virtual void initLayout();

    ///Read settings
  void readSettings();
  
  /**@name Utility functions */
  //@{
  /// Load the user file specified in the text field
  bool loadUserFile();

  /// Read a limit line from the user file
  void readLimits(const QString & com_line);
  //@}

  /**@name Python code utility commands */
  //@{
  /// Construct a LoadRaw Python command
  QString writeLoadRawCmd(const QString & filename, const QString & workspace, const QString & spec_min = QString(), 
			  const QString & spec_max = QString(), const QString & spec_list = QString(), 
			  const QString & cache_opt = QString());
  //@}
			   
private slots:
  /// Select the data directory
  void selectDataDir();

  /// Select the user file
  void selectUserFile();

  /// Receive a load button click
  void loadButtonClicked();

  /// Plot button has been clicked
  void plotButtonClicked();

  /// A ComboBox option change
  void stepComboChange(int new_index);

  ///Save settings
  void saveSettings();

private:
  /// The form generated by Qt Designer
  Ui::SANSRunWindow m_uiForm;

  /// The last directory that was viewed
  QString m_last_dir;
  
  /// A map for quickly retrieving the different line edits
  QMap<int, QLineEdit*> m_run_no_boxes;

  /// A list of unique run numbers that have been loaded
  QStringList m_unique_runs;

  /// A hash for quickly retrieving the different label fields
  QHash<int, QLabel*> m_period_lbls;
};

}
}

#endif //MANTIDQTCUSTOMINTERFACES_SANSRUNWINDOW_H_
