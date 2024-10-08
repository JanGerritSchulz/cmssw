#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/DDCMS/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/plugins/CmsDetConstruction.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <memory>

template <class FilteredView>
void CmsDetConstruction<FilteredView>::buildSmallDetsforGlued(FilteredView& fv,
                                                              GeometricDet* mother,
                                                              const std::string& attribute) {
  auto det = std::make_unique<GeometricDet>(&fv,
                                            CmsTrackerLevelBuilder<FilteredView>::theCmsTrackerStringToEnum.type(
                                                ExtractStringFromDDD<FilteredView>::getString(attribute, &fv)));
  if (det->stereo()) {
    uint32_t temp = 1;
    det->setGeographicalID(DetId(temp));
  } else {
    uint32_t temp = 2;
    det->setGeographicalID(DetId(temp));
  }

  mother->addComponent(det.release());
}

template <class FilteredView>
void CmsDetConstruction<FilteredView>::buildSmallDetsforStack(FilteredView& fv,
                                                              GeometricDet* mother,
                                                              const std::string& attribute) {
  auto det = std::make_unique<GeometricDet>(&fv,
                                            CmsTrackerLevelBuilder<FilteredView>::theCmsTrackerStringToEnum.type(
                                                ExtractStringFromDDD<FilteredView>::getString(attribute, &fv)));

  if (det->isLowerSensor()) {
    uint32_t temp = 1;
    det->setGeographicalID(DetId(temp));
  } else if (det->isUpperSensor()) {
    uint32_t temp = 2;
    det->setGeographicalID(DetId(temp));
  } else {
    edm::LogError("DetConstruction") << " module defined in a Stack but not upper either lower!? ";
  }
  mother->addComponent(det.release());
}

template <class FilteredView>
void CmsDetConstruction<FilteredView>::buildSmallDetsfor3D(FilteredView& fv,
                                                           GeometricDet* mother,
                                                           const std::string& attribute) {
  auto det = std::make_unique<GeometricDet>(&fv,
                                            CmsTrackerLevelBuilder<FilteredView>::theCmsTrackerStringToEnum.type(
                                                ExtractStringFromDDD<FilteredView>::getString(attribute, &fv)));

  if (det->isFirstSensor()) {
    uint32_t temp = 1;
    det->setGeographicalID(DetId(temp));
  } else if (det->isSecondSensor()) {
    uint32_t temp = 2;
    det->setGeographicalID(DetId(temp));
  } else {
    edm::LogError("DetConstruction") << " module defined in a 3D module but not first or second sensor!? ";
  }
  mother->addComponent(det.release());
}

/*
 * OLD DD.
 * Module with 2 sensors: calculate the sensor local ID, and add the sensor to its mother volume (module).
 * Module with 1 sensor: just add the sensor to to its mother volume (ladder).
 */
template <>
void CmsDetConstruction<DDFilteredView>::buildComponent(DDFilteredView& fv,
                                                        GeometricDet* mother,
                                                        const std::string& attribute) {
  // Mother volume
  // Module with 2 sensors: the mother volume is the module volume.
  // Module with 1 sensor: the mother volume is the ladder volume.
  const std::string& myTopologicalNameInXMLs = ExtractStringFromDDD<DDFilteredView>::getString(attribute, &fv);
  const GeometricDet::GDEnumType& myTopologicalType =
      CmsTrackerLevelBuilder<DDFilteredView>::theCmsTrackerStringToEnum.type(myTopologicalNameInXMLs);

  auto det = std::make_unique<GeometricDet>(&fv, myTopologicalType);

  const bool isPhase1ModuleWith2Sensors = (myTopologicalType == GeometricDet::mergedDet);
  const bool isPhase2ModuleWith2Sensors = (myTopologicalType == GeometricDet::OTPhase2Stack);
  const bool isPhase2BarrelModuleWith2Sensors = (myTopologicalType == GeometricDet::ITPhase2Combined);

  // CASE A: MODULE HAS 2 SENSORS
  if (isPhase1ModuleWith2Sensors || isPhase2ModuleWith2Sensors || isPhase2BarrelModuleWith2Sensors) {
    // Go down in hierarchy: from module to sensor
    bool dodets = fv.firstChild();  // very important
    while (dodets) {
      // PHASE 1 (MERGEDDET)
      if (isPhase1ModuleWith2Sensors) {
        buildSmallDetsforGlued(fv, det.get(), attribute);
      }
      // PHASE 2 (STACKDET)
      else if (isPhase2ModuleWith2Sensors) {
        buildSmallDetsforStack(fv, det.get(), attribute);
      } else if (isPhase2BarrelModuleWith2Sensors) {
        buildSmallDetsfor3D(fv, det.get(), attribute);
      }

      dodets = fv.nextSibling();
    }

    fv.parent();
  }

  // CASE B: MODULE HAS 1 SENSOR: NOTHING SPECIFIC TO DO
  // Indeed, we are not going to sort sensors within module, if there is only 1 sensor!

  // ALL CASES: add sensor to its mother volume (module or ladder).
  mother->addComponent(det.release());
}

/*
 * DD4hep.
 * Module with 2 sensors: calculate the sensor local ID, and add the sensor to its mother volume (module).
 * Module with 1 sensor: just add the sensor to its mother volume (ladder).
 */
template <>
void CmsDetConstruction<cms::DDFilteredView>::buildComponent(cms::DDFilteredView& fv,
                                                             GeometricDet* mother,
                                                             const std::string& attribute) {
  // Mother volume
  // Module with 2 sensors: the mother volume is the module volume.
  // Module with 1 sensor: the mother volume is the ladder volume.
  const std::string& myTopologicalNameInXMLs = ExtractStringFromDDD<cms::DDFilteredView>::getString(attribute, &fv);
  const GeometricDet::GDEnumType& myTopologicalType =
      CmsTrackerLevelBuilder<cms::DDFilteredView>::theCmsTrackerStringToEnum.type(myTopologicalNameInXMLs);
  auto det = std::make_unique<GeometricDet>(&fv, myTopologicalType);

  const bool isPhase1ModuleWith2Sensors = (myTopologicalType == GeometricDet::mergedDet);
  const bool isPhase2ModuleWith2Sensors = (myTopologicalType == GeometricDet::OTPhase2Stack);
  const bool isPhase2BarrelModuleWith2Sensors = (myTopologicalType == GeometricDet::ITPhase2Combined);

  // CASE A: MODULE HAS 2 SENSORS
  if (isPhase1ModuleWith2Sensors || isPhase2ModuleWith2Sensors || isPhase2BarrelModuleWith2Sensors) {
    // Go down in hierarchy: from module to sensor
    if (!fv.firstChild()) {  // very important
      edm::LogError("CmsDetConstruction::buildComponent. Cannot go down to sensor volume.");
      return;
    }

    // This is the sensor hierarchy level
    const int sensorHierarchyLevel = fv.level();

    // Loop on all siblings (ie, on all sensors)
    while (fv.level() == sensorHierarchyLevel) {
      // PHASE 1 (MERGEDDET)
      if (isPhase1ModuleWith2Sensors) {
        buildSmallDetsforGlued(fv, det.get(), attribute);
      }
      // PHASE 2 (STACKDET)
      else if (isPhase2ModuleWith2Sensors) {
        buildSmallDetsforStack(fv, det.get(), attribute);
      } else if (isPhase2BarrelModuleWith2Sensors) {
        buildSmallDetsfor3D(fv, det.get(), attribute);
      }

      // Go to the next volume in FilteredView.
      // NB: If this volume is another sensor of the same module, will stay in the loop.
      // Otherwise, it is very important to access the next volume to be treated anyway.
      fv.firstChild();
    }
  }

  // CASE B: MODULE HAS 1 SENSOR: NOTHING SPECIFIC TO DO
  // Indeed, we are not going to sort sensors within module, if there is only 1 sensor!
  else {
    // Go to the next volume in FilteredView.
    fv.firstChild();
  }

  // ALL CASES: add sensor to its mother volume (module or ladder).
  mother->addComponent(det.release());
}
