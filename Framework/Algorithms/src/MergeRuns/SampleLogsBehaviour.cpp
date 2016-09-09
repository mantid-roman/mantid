#include <MantidAlgorithms/MergeRuns/SampleLogsBehaviour.h>
#include "MantidGeometry/Instrument.h"
#include "MantidKernel/TimeSeriesProperty.h"

namespace Mantid {
namespace Algorithms {

using namespace Kernel;
using namespace API;

namespace {
std::string generateDifferenceMessage(const std::string &item,
                                      const std::string &wsName,
                                      const std::string &wsValue,
                                      const std::string &firstValue) {
  std::stringstream stringstream;
  stringstream << "Item \"" << item
               << "\" has different values in workspaces! Found: " << wsValue
               << " in workspace " << wsName
               << " but value in first workspace value was: " << firstValue
               << "." << std::endl;
  return stringstream.str();
}
}

const std::string SampleLogsBehaviour::TIME_SERIES_MERGE =
    "sample_logs_time_series";
const std::string SampleLogsBehaviour::LIST_MERGE = "sample_logs_list";
const std::string SampleLogsBehaviour::WARN_MERGE = "sample_logs_warn";
const std::string SampleLogsBehaviour::FAIL_MERGE = "sample_logs_fail";
const std::string SampleLogsBehaviour::WARN_MERGE_TOLERANCES =
    "sample_logs_warn_tolerances";
const std::string SampleLogsBehaviour::FAIL_MERGE_TOLERANCES =
    "sample_logs_fail_tolerances";

const std::string SampleLogsBehaviour::TIME_SERIES_SUFFIX = "_time_series";
const std::string SampleLogsBehaviour::LIST_SUFFIX = "_list";

SampleLogsBehaviour::SampleLogsBehaviour(
    const MatrixWorkspace_sptr &ws, Logger &g_log,
    std::string sampleLogsTimeSeries, std::string sampleLogsList,
    std::string sampleLogsWarn, std::string sampleLogsWarnTolerances,
    std::string sampleLogsFail, std::string sampleLogsFailTolerances)
    : m_logger(g_log) {
  updateSampleMap(m_logMap, time_series, sampleLogsTimeSeries, ws, "");
  updateSampleMap(m_logMap, list, sampleLogsList, ws, "");
  updateSampleMap(m_logMap, warn, sampleLogsWarn, ws, sampleLogsWarnTolerances);
  updateSampleMap(m_logMap, fail, sampleLogsFail, ws, sampleLogsFailTolerances);

  SampleLogsMap instrumentMap;
  this->createSampleLogsMapsFromInstrumentParams(instrumentMap, ws);

  // This adds the parameters from the instrument to the main map, with any
  // duplicates left as the versions in the MergeRuns arguments.
  m_logMap.insert(instrumentMap.begin(), instrumentMap.end());
}

void SampleLogsBehaviour::createSampleLogsMapsFromInstrumentParams(
    SampleLogsMap &map, const MatrixWorkspace_sptr &ws) {
  std::string params =
      ws->getInstrument()->getParameterAsString(TIME_SERIES_MERGE, false);
  updateSampleMap(map, time_series, params, ws, "", true);

  params = ws->getInstrument()->getParameterAsString(LIST_MERGE, false);
  updateSampleMap(map, list, params, ws, "", true);

  params = ws->getInstrument()->getParameterAsString(WARN_MERGE, false);
  std::string paramsTolerances;
  paramsTolerances =
      ws->getInstrument()->getParameterAsString(WARN_MERGE_TOLERANCES, false);
  updateSampleMap(map, warn, params, ws, paramsTolerances, true);

  params = ws->getInstrument()->getParameterAsString(FAIL_MERGE, false);
  paramsTolerances =
      ws->getInstrument()->getParameterAsString(FAIL_MERGE_TOLERANCES, false);
  updateSampleMap(map, fail, params, ws, paramsTolerances, true);
}

void SampleLogsBehaviour::updateSampleMap(
    SampleLogsMap &map, const MergeLogType &sampleLogBehaviour,
    const std::string &params, const API::MatrixWorkspace_sptr &ws,
    const std::string paramsTolerances, bool skipIfInPrimaryMap) {

  StringTokenizer tokenizer(params, ",", StringTokenizer::TOK_TRIM |
                                             StringTokenizer::TOK_IGNORE_EMPTY);
  StringTokenizer tokenizerTolerances(paramsTolerances, ",",
                                      StringTokenizer::TOK_TRIM |
                                          StringTokenizer::TOK_IGNORE_EMPTY);

  auto tolerancesStringVector = tokenizerTolerances.asVector();

  std::vector<double> tolerancesVector = createTolerancesVector(
      tokenizer.asVector(), tolerancesStringVector, paramsTolerances);

  StringTokenizer::Iterator i = tokenizer.begin();
  std::vector<double>::iterator j = tolerancesVector.begin();

  for (; i != tokenizer.end() && j != tolerancesVector.end(); ++i, ++j) {
    auto item = *i;
    auto tolerance = *j;

    if (skipIfInPrimaryMap && (m_logMap.count(item) != 0)) {
      continue;
    }

    if (map.count(item) != 0) {
      throw std::invalid_argument(
          "Error when making list of merge items, sample log \"" + item +
          "\" defined more than once!");
    }

    std::shared_ptr<Property> prop;
    try {
      prop = std::shared_ptr<Property>(ws->getLog(item)->clone());
    } catch (std::invalid_argument &) {
      m_logger.warning() << "Could not merge sample log \"" << item
                         << "\", does not exist in workspace!" << std::endl;
      continue;
    }

    bool isNumeric;
    double value = 0.0;
    try {
      value = ws->getLogAsSingleValue(item);
      isNumeric = true;
    } catch (std::invalid_argument &) {
      isNumeric = false;
      if (sampleLogBehaviour == time_series) {
        m_logger.error() << item << " could not be converted to a numeric type";
        continue;
      }
    }

    if (sampleLogBehaviour == time_series) {
      try {
        // See if property exists already - merging an output of MergeRuns
        prop = std::shared_ptr<Property>(
            ws->getLog(item + TIME_SERIES_SUFFIX)->clone());
      } catch (std::invalid_argument &) {
        // Property does not already exist
        std::unique_ptr<Kernel::TimeSeriesProperty<double>> timeSeriesProp(
            new TimeSeriesProperty<double>(item + TIME_SERIES_SUFFIX));
        std::string startTime = ws->mutableRun().startTime().toISO8601String();
        timeSeriesProp->addValue(startTime, value);
        ws->mutableRun().addLogData(
            std::unique_ptr<Kernel::Property>(std::move(timeSeriesProp)));
        prop = std::shared_ptr<Property>(
            ws->getLog(item + TIME_SERIES_SUFFIX)->clone());
      }
    } else if (sampleLogBehaviour == list) {
      try {
        // See if property exists already - merging an output of MergeRuns
        prop =
            std::shared_ptr<Property>(ws->getLog(item + LIST_SUFFIX)->clone());
      } catch (std::invalid_argument &) {
        ws->mutableRun().addProperty(item + LIST_SUFFIX, prop->value());
        prop =
            std::shared_ptr<Property>(ws->getLog(item + LIST_SUFFIX)->clone());
      }
    }

    map[item] = {sampleLogBehaviour, prop, tolerance, isNumeric};
  }
}

std::vector<double> SampleLogsBehaviour::createTolerancesVector(
    const std::vector<std::string> &names,
    const std::vector<std::string> &tolerances,
    const std::string paramsTolerances) {
  size_t numberNames = names.size();
  size_t numberTolerances = tolerances.size();

  std::vector<double> tolerancesVector(numberNames);

  if (numberNames == numberTolerances) {
    std::transform(tolerances.begin(), tolerances.end(),
                   tolerancesVector.begin(),
                   [](const std::string &val) { return std::stod(val); });
  } else if (paramsTolerances.empty()) {
    std::fill(tolerancesVector.begin(), tolerancesVector.end(), -1.0);
  } else if (numberTolerances == 1) {
    double value = std::stod(tolerances.front());
    std::fill(tolerancesVector.begin(), tolerancesVector.end(), value);
  } else {
    throw std::invalid_argument("Invalid length of tolerances, found " +
                                std::to_string(numberTolerances) +
                                " tolerance values but " +
                                std::to_string(numberNames) + " names.");
  }

  return tolerancesVector;
}

void SampleLogsBehaviour::calculateUpdatedSampleLogs(
    const MatrixWorkspace_sptr &addeeWS,
    const API::MatrixWorkspace_sptr &outWS) {
  for (auto item : m_logMap) {
    Property *addeeWSProperty = addeeWS->getLog(item.first);

    double addeeWSNumber = 0;
    double outWSNumber = 0;

    try {
      addeeWSNumber = addeeWS->getLogAsSingleValue(item.first);
      outWSNumber = outWS->getLogAsSingleValue(item.first);
    } catch (std::invalid_argument &) {
      if (item.second.isNumeric) {
        throw std::invalid_argument(
            item.first + " could not be converted to a numeric type");
      }
    }

    switch (item.second.type) {
    case time_series: {
      this->updateTimeSeriesProperty(addeeWS, outWS, item.first);
      break;
    }
    case list: {
      this->updateListProperty(addeeWS, outWS, addeeWSProperty, item.first);
      break;
    }
    case warn:
      this->checkWarnProperty(addeeWS, addeeWSProperty, item.second,
                              addeeWSNumber, outWSNumber, item.first);
      break;
    case fail:
      this->checkErrorProperty(addeeWS, addeeWSProperty, item.second,
                               addeeWSNumber, outWSNumber, item.first);
      break;
    }
  }
}

void SampleLogsBehaviour::updateTimeSeriesProperty(
    const MatrixWorkspace_sptr &addeeWS, const MatrixWorkspace_sptr &outWS,
    const std::string name) {
  try {
    // If this already exists we do not need to do anything, Time Series Logs
    // are combined when adding workspaces.
    addeeWS->getLog(name + TIME_SERIES_SUFFIX);
  } catch (std::invalid_argument &) {
    auto timeSeriesProp = outWS->mutableRun().getTimeSeriesProperty<double>(
        name + TIME_SERIES_SUFFIX);
    Kernel::DateAndTime startTime = addeeWS->mutableRun().startTime();
    double value = addeeWS->mutableRun().getLogAsSingleValue(name);
    timeSeriesProp->addValue(startTime, value);
  }
}

void SampleLogsBehaviour::updateListProperty(
    const MatrixWorkspace_sptr &addeeWS, const MatrixWorkspace_sptr &outWS,
    Property *addeeWSProperty, const std::string name) {
  try {
    // If this already exists we combine the two strings.
    auto propertyAddeeWS = addeeWS->getLog(name + LIST_SUFFIX);
    auto propertyOutWS = outWS->mutableRun().getProperty(name + LIST_SUFFIX);
    propertyOutWS->setValue(propertyOutWS->value() + ", " +
                            propertyAddeeWS->value());
  } catch (std::invalid_argument &) {
    auto property = outWS->mutableRun().getProperty(name + LIST_SUFFIX);
    property->setValue(property->value() + ", " + addeeWSProperty->value());
  }
}

void SampleLogsBehaviour::checkWarnProperty(const MatrixWorkspace_sptr &addeeWS,
                                            Property *addeeWSProperty,
                                            const SampleLogBehaviour &behaviour,
                                            const double addeeWSNumber,
                                            const double outWSNumber,
                                            const std::string name) {
  if (behaviour.isNumeric && behaviour.tolerance > 0.0) {
    if (std::abs(addeeWSNumber - outWSNumber) > behaviour.tolerance) {
      m_logger.warning() << generateDifferenceMessage(
          name, addeeWS->name(), addeeWSProperty->value(),
          behaviour.property->value());
    }
  } else {
    if (behaviour.property->value().compare(addeeWSProperty->value()) != 0) {
      m_logger.warning() << generateDifferenceMessage(
          name, addeeWS->name(), addeeWSProperty->value(),
          behaviour.property->value());
    }
  }
}

void SampleLogsBehaviour::checkErrorProperty(
    const MatrixWorkspace_sptr &addeeWS, Property *addeeWSProperty,
    const SampleLogBehaviour &behaviour, const double addeeWSNumber,
    const double outWSNumber, const std::string name) {
  if (behaviour.isNumeric && behaviour.tolerance > 0.0) {
    if (std::abs(addeeWSNumber - outWSNumber) > behaviour.tolerance) {
      throw std::invalid_argument(generateDifferenceMessage(
          name, addeeWS->name(), addeeWSProperty->value(),
          behaviour.property->value()));
    }
  } else {
    if (behaviour.property->value().compare(addeeWSProperty->value()) != 0) {
      throw std::invalid_argument(generateDifferenceMessage(
          name, addeeWS->name(), addeeWSProperty->value(),
          behaviour.property->value()));
    }
  }
}

void SampleLogsBehaviour::setUpdatedSampleLogs(
    const API::MatrixWorkspace_sptr &ws) {
  for (auto item : m_logMap) {
    std::string propertyToReset = item.first;

    if (item.second.type == time_series) {
      propertyToReset = item.first + TIME_SERIES_SUFFIX;
    } else if (item.second.type == list) {
      propertyToReset = item.first + LIST_SUFFIX;
    } else {
      return;
    }

    Property *outWSProperty = ws->mutableRun().getProperty(propertyToReset);
    item.second.property = std::shared_ptr<Property>(outWSProperty->clone());
  }
}

void SampleLogsBehaviour::resetSampleLogs(const API::MatrixWorkspace_sptr &ws) {
  for (auto item : m_logMap) {
    std::string propertyToReset = item.first;

    if (item.second.type == time_series) {
      propertyToReset = item.first + TIME_SERIES_SUFFIX;
      auto property =
          std::unique_ptr<Kernel::Property>(item.second.property->clone());
      ws->mutableRun().addProperty(std::move(property), true);
    } else if (item.second.type == list) {
      propertyToReset = item.first + LIST_SUFFIX;
      ws->mutableRun()
          .getProperty(propertyToReset)
          ->setValue(item.second.property->value());
    } else {
      return;
    }
  }
}

} // namespace Algorithms
} // namespace Mantid
