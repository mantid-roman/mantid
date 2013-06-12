/*WIKI* 

Sets the neutrons information in the sample. You can either enter details about the chemical formula or atomic number, 
or you can provide specific values for the attenuation and scattering cross sections and the sample number density.  
If you decide to provide specific values you must give values for all three (attenuation and scattering cross sections and the sample number density), and any formula information will be ignored.
If you miss any of the three specific values then the other will be ignored.

Neutron scattering lengths and cross sections of the elements and their isotopes have been taken from [http://www.ncnr.nist.gov/resources/n-lengths/list.html].
 *WIKI*/
/*WIKI_USAGE* 
=====Setting the sample by simple formula=====
 SetSampleMaterial(InputWorkspace='IRS26173',ChemicalFormula='Fe')
 
=====Setting the sample by a more complex formula=====
 SetSampleMaterial(InputWorkspace='IRS26173',ChemicalFormula='Al2-O3', UnitCellVolume='253.54', ZParameter='6')

=====Setting the sample by specific values (all three must be specified)=====
 SetSampleMaterial(InputWorkspace='IRS26173',AttenuationXSection=2.56,ScatteringXSection=11.62,SampleNumberDensity=0.0849106)

=====Extracting the set values out by python=====
 sam = ws.sample()
 mat = sam.getMaterial()
 print mat.absorbXSection()
  1.3374
 print mat.cohScatterXSection()
  339.1712
 print mat.name()
  C2 H4
 print mat.totalScatterXSection()
  339.1712

 *WIKI_USAGE*/
//--------------------------------
// Includes
//--------------------------------
#include "MantidDataHandling/SetSampleMaterial.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/Sample.h"
#include "MantidGeometry/Crystal/OrientedLattice.h"
#include "MantidKernel/MandatoryValidator.h"
#include "MantidKernel/Atom.h"
#include "MantidKernel/EnabledWhenProperty.h"
#include "MantidKernel/NeutronAtom.h"
#include "MantidKernel/Material.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/PhysicalConstants.h"

using namespace Mantid::PhysicalConstants;

namespace Mantid
{
namespace DataHandling
{
  // Register the algorithm into the AlgorithmFactory
  DECLARE_ALGORITHM(SetSampleMaterial)

  /// Sets documentation strings for this algorithm
  void SetSampleMaterial::initDocs()
  {
    this->setWikiSummary("Sets the neutrons information in the sample.");
    this->setOptionalMessage("Sets the neutrons information in the sample.");
  }

  using namespace Mantid::DataHandling;
  using namespace Mantid::API;
  using namespace Kernel;

  void SetSampleMaterial::logMaterial(const Material *mat)
  {
    g_log.notice() << "Sample number density ";
    if (mat->numberDensity() == 1.)
      g_log.notice() << "(NOT SPECIFIED) ";
    g_log.notice() << "= " << mat->numberDensity() << " atoms/Angstrom^3\n";
    g_log.notice() << "Cross sections for wavelength = " << NeutronAtom::ReferenceLambda << "Angstroms\n"
                   << "Coherent "   << mat->cohScatterXSection() << " barns\n"
                   << "Incoherent " << mat->incohScatterXSection() << " barns\n"
                   << "Total "      << mat->totalScatterXSection() << " barns\n"
                   << "Absorption " << mat->absorbXSection() << " barns\n";
  }

  /**
   * Initialize the algorithm
   */
  void SetSampleMaterial::init()
  {
    using namespace Mantid::Kernel;
    declareProperty(
        new WorkspaceProperty<MatrixWorkspace>("InputWorkspace","",Direction::InOut),
        "The workspace with which to associate the sample ");
    declareProperty("ChemicalFormula", "", "ChemicalFormula or AtomicNumber must be given. "
        "Enter a composition as a molecular formula of \n"
        "elements or isotopes.  For example, basic "
        "elements might be \"H\", \"Fe\" or \"Si\", etc. \n"
        "A molecular formula of elements might be "
        "\"H4-N2-C3\", which corresponds to a molecule \n"
        "with 4 Hydrogen atoms, 2 Nitrogen atoms and "
        "3 Carbon atoms.  Each element in a molecular \n"
        "formula is followed by the number of the atoms "
        "for that element, specified _without a hyphen_, \n"
        "because each element is separated from other "
        "elements using a hyphen.  The number of atoms \n"
        "can be integer or float, but must start with "
        "a digit, e.g. 0.6 is fine but .6 is not. \n"
        "Isotopes may also be included in a material "
        "composition, and can be specified alone \n"
        "(as in \"Li7\"), or in a molecular formula "
        "(as in \"(Li7)2-C-H4-N-Cl6\").  Note, however, \n"
        "that No Spaces or Hyphens are allowed in an "
        "isotope symbol specification.  Also Note \n"
        "that for isotopes specified in a molecular "
        "expression, the isotope must be enclosed \n"
        "by parenthesis, except for two special "
        "cases, D and T, which stand for H2 and H3, \n"
        "respectively.");
    declareProperty("AtomicNumber", 0, "ChemicalFormula or AtomicNumber must be given");
    declareProperty("MassNumber", 0, "Mass number if ion (default is 0)");
    auto mustBePositive = boost::make_shared<BoundedValidator<double> >();
    mustBePositive->setLower(0.0);
    declareProperty("SampleNumberDensity", EMPTY_DBL(), mustBePositive,
        "Optional:  This number density of the sample in number of atoms per cubic angstrom will be used instead of calculated");
    declareProperty("ZParameter", EMPTY_DBL(), mustBePositive,
        "Number of formulas in the unit cell needed for chemical formulas with more than 1 atom");
    declareProperty("UnitCellVolume", EMPTY_DBL(), mustBePositive,
        "Unit cell volume in Angstoms^3 needed for chemical formulas with more than 1 atom");
    declareProperty("AttenuationXSection", EMPTY_DBL(), mustBePositive,
        "Optional:  This absorption cross-section for the sample material in barns will be used instead of calculated");
    declareProperty("ScatteringXSection", EMPTY_DBL(), mustBePositive,
        "Optional:  This scattering cross-section (coherent + incoherent) for the sample material in barns will be used instead of calculated");

	
    // Perform Group Associations.
    std::string formulaGrp("By Formula or Atomic Number");
    setPropertyGroup("ChemicalFormula", formulaGrp);
    setPropertyGroup("AtomicNumber", formulaGrp);
    setPropertyGroup("MassNumber", formulaGrp);

    std::string densityGrp("Sample Density");
    setPropertyGroup("SampleNumberDensity", densityGrp);
    setPropertyGroup("ZParameter", densityGrp);
    setPropertyGroup("UnitCellVolume", densityGrp);

    std::string specificValuesGrp("Enter Specific Values");
    setPropertyGroup("AttenuationXSection", specificValuesGrp);
    setPropertyGroup("ScatteringXSection", specificValuesGrp);

    // Extra property settings
    setPropertySettings("AtomicNumber", new Kernel::EnabledWhenProperty("ChemicalFormula", Kernel::IS_DEFAULT));
    setPropertySettings("MassNumber", new Kernel::EnabledWhenProperty("ChemicalFormula", Kernel::IS_DEFAULT));
    setPropertySettings("UnitCellVolume", new Kernel::EnabledWhenProperty("SampleNumberDensity", Kernel::IS_DEFAULT));
    setPropertySettings("ZParameter", new Kernel::EnabledWhenProperty("SampleNumberDensity", Kernel::IS_DEFAULT));

  }

  /**
   * Execute the algorithm
   */
  void SetSampleMaterial::exec()
  {
    // Get the input workspace
    MatrixWorkspace_sptr workspace = getProperty("InputWorkspace");

    // determine the sample number density
    double rho = getProperty("SampleNumberDensity"); // in Angstroms-3
    if (isEmpty(rho))
    {
      // if rho isn't set then just choose it to be one
      rho = 1.;

      double unitCellVolume = getProperty("UnitCellVolume"); // in Angstroms^3
      double zParameter = getProperty("ZParameter"); // number of atoms

      // get the unit cell volume from the workspace if it isn't set
      if (isEmpty(unitCellVolume) && workspace->sample().hasOrientedLattice())
      {
        unitCellVolume = workspace->sample().getOrientedLattice().volume();
        g_log.notice() << "found unit cell volume " << unitCellVolume << " Angstrom^-3\n";
      }
      // density is just number of atoms in the unit cell
      // ...but only calculate it if you have both numbers
      if ((!isEmpty(zParameter)) && (!isEmpty(unitCellVolume)))
        rho = zParameter / unitCellVolume;
    }

    // get the scattering information
    const std::string chemicalSymbol = getProperty("ChemicalFormula");
    const int z_number = getProperty("AtomicNumber");
    const int a_number = getProperty("MassNumber");
    double sigma_atten = getProperty("AttenuationXSection"); // in barns
    double sigma_s = getProperty("ScatteringXSection"); // in barns

    // Use user variables if all three are given
    if (sigma_atten != EMPTY_DBL() && sigma_s != EMPTY_DBL() && rho != EMPTY_DBL())
    {
      NeutronAtom *neutron = new NeutronAtom(static_cast<uint16_t>(z_number), static_cast<uint16_t>(a_number),
                                             0.0, 0.0, sigma_s, 0.0, sigma_s, sigma_atten);
      Material *mat = new Material(chemicalSymbol, *neutron, rho);
      workspace->mutableSample().setMaterial(*mat);
      logMaterial(mat);
      return;
    }

    // Use chemical symbol if given by user
    try
    {
      Atom myAtom = getAtom(chemicalSymbol, static_cast<uint16_t>(a_number));
      Material *mat = new Material(chemicalSymbol, myAtom.neutron, myAtom.number_density);
      workspace->mutableSample().setMaterial(*mat);
      logMaterial(mat);
    }
    catch (...)
    {
      // Use chemical formula if given by user
      try
      {
        Material::ChemicalFormula CF = Material::parseChemicalFormula(chemicalSymbol);
        g_log.notice() << "Found " << CF.atoms.size() << " atoms in \"" << chemicalSymbol << "\"\n";

        double numAtoms = 0.; // number of atoms in formula
        NeutronAtom neutron(0, 0., 0., 0., 0., 0., 0.); // starting thing for neutronic information
        for (size_t i=0; i<CF.atoms.size(); i++)
        {
          Atom myAtom = getAtom(CF.atoms[i], CF.aNumbers[i]);
          neutron = neutron + CF.numberAtoms[i] * myAtom.neutron;

          g_log.information() << myAtom << ": " << myAtom.neutron << "\n";
          numAtoms += static_cast<double>(CF.numberAtoms[i]);
        }
        // normalize the accumulated number by the number of atoms
        neutron = (1. / numAtoms) * neutron; // funny syntax b/c of operators in neutron atom

        // create the material
        Material *mat = new Material(chemicalSymbol, neutron, rho);
        workspace->mutableSample().setMaterial(*mat);
        logMaterial(mat);
      }
      catch (...)
      {
        // Use atomic and mass number if chemical formula does not work
        try
        {
          Atom myAtom = getAtom(static_cast<uint16_t>(z_number), static_cast<uint16_t>(a_number));
          Material *mat = new Material(chemicalSymbol, myAtom.neutron, myAtom.number_density);
          workspace->mutableSample().setMaterial(*mat);
          logMaterial(mat);
        }
        catch(std::invalid_argument&)
        {
          g_log.notice("ChemicalFormula or AtomicNumber was not found in table.");
          throw std::invalid_argument("ChemicalFormula or AtomicNumber was not found in table");
        }
      }
    }
    // Done!
    progress(1);
  }

}
}
