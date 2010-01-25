#ifndef PROPERTYMANAGERTEST_H_
#define PROPERTYMANAGERTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidKernel/PropertyManager.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/MandatoryValidator.h"

using namespace Mantid::Kernel;

class PropertyManagerHelper : public PropertyManager
{
public:
  PropertyManagerHelper() : PropertyManager() {}

  using PropertyManager::declareProperty;
  using PropertyManager::setProperty;
  using PropertyManager::getPointerToProperty;
};

class PropertyManagerTest : public CxxTest::TestSuite
{
public:
  PropertyManagerTest()
  {
    Property *p = new PropertyWithValue<int>("aProp",1);
    manager.declareProperty(p);
    manager.declareProperty("anotherProp",1.11);
    manager.declareProperty("yetAnotherProp","itsValue");
  }

  void testConstructor()
  {
    PropertyManagerHelper mgr;
    std::vector<Property*> props = mgr.getProperties();
    TS_ASSERT( props.empty() );
  }

  void testdeclareProperty_pointer()
  {
    PropertyManagerHelper mgr;
    Property *p = new PropertyWithValue<double>("myProp", 9.99);
    TS_ASSERT_THROWS_NOTHING( mgr.declareProperty(p) );
    TS_ASSERT( mgr.existsProperty(p->name()) );
    // Confirm that the first 4 characters of the string are the same
    
    // Note that some versions of boost::lexical_cast > 1.34 give a string such as
    // 9.9900000000000002 rather than 9.99. Converting back to a double however does
    // still give the correct 9.99.
    
    TS_ASSERT_EQUALS( mgr.getPropertyValue("myProp").substr(0,4),  std::string("9.99") );
      
    TS_ASSERT_THROWS( mgr.declareProperty(p), Exception::ExistsError );
    TS_ASSERT_THROWS( mgr.declareProperty(new PropertyWithValue<int>("",0)), std::invalid_argument );
    mgr.declareProperty(new PropertyWithValue<int>("GoodIntProp",1), "Test doc"); 
    TS_ASSERT_EQUALS( mgr.getPointerToProperty("GoodIntProp")->documentation(), "Test doc" );
  }

  void testdeclareProperty_int()
  {
    PropertyManagerHelper mgr;
    TS_ASSERT_THROWS_NOTHING( mgr.declareProperty("myProp", 1) )
    TS_ASSERT( ! mgr.getPropertyValue("myProp").compare("1") )
    TS_ASSERT_THROWS( mgr.declareProperty("MYPROP", 5), Exception::ExistsError )
    TS_ASSERT_THROWS( mgr.declareProperty("", 5), std::invalid_argument )
  }

  void testdeclareProperty_double()
  {
    PropertyManagerHelper mgr;
    BoundedValidator<double> *v = new BoundedValidator<double>(1,5);
    TS_ASSERT_THROWS_NOTHING( mgr.declareProperty("myProp", 9.99, v) )
    // Note that some versions of boost::lexical_cast > 1.34 give a string such as
    // 9.9900000000000002 rather than 9.99. Converting back to a double however does
    // still give the correct 9.99.
    TS_ASSERT_EQUALS( mgr.getPropertyValue("myProp").substr(0,4), std::string("9.99") )
    TS_ASSERT_THROWS_NOTHING( mgr.declareProperty("withDoc", 4.4, v->clone(), "Test doc doub") )
    TS_ASSERT_EQUALS( mgr.getPointerToProperty("withDoc")->documentation(), "Test doc doub" )
    TS_ASSERT_THROWS( mgr.declareProperty("MYPROP", 5.5), Exception::ExistsError )
    TS_ASSERT_THROWS( mgr.declareProperty("", 5.5), std::invalid_argument )
  }

  void testdeclareProperty_string()
  {
    PropertyManagerHelper mgr;
    TS_ASSERT_THROWS_NOTHING( mgr.declareProperty("myProp", "theValue", new MandatoryValidator<std::string>, "hello") )
    TS_ASSERT_EQUALS( mgr.getPropertyValue("myProp"), "theValue" )
      Property *p;
    TS_ASSERT_THROWS_NOTHING( p = mgr.getProperty("myProp") )
    TS_ASSERT_EQUALS(p->documentation(),"hello");

    TS_ASSERT_THROWS( mgr.declareProperty("MYPROP", "aValue"), Exception::ExistsError )
    TS_ASSERT_THROWS( mgr.declareProperty("", "aValue"), std::invalid_argument )
  }

  void testSetProperties()
  {
    PropertyManagerHelper mgr;
    mgr.declareProperty("APROP", 1);
    mgr.declareProperty("anotherProp", 1.0);
    TS_ASSERT_THROWS_NOTHING( mgr.setProperties("APROP=15;anotherProp=1.3") )
    TS_ASSERT( ! mgr.getPropertyValue("APROP").compare("15") )
    TS_ASSERT( ! mgr.getPropertyValue("anotherProp").compare("1.3") )
  }

  void testSetPropertyValue()
  {
    manager.setPropertyValue("APROP","10");
    TS_ASSERT( ! manager.getPropertyValue("aProp").compare("10") )
    manager.setPropertyValue("aProp","1");
    TS_ASSERT_THROWS( manager.setPropertyValue("fhfjsdf","0"), Exception::NotFoundError )
  }

  void testSetProperty()
  {
    TS_ASSERT_THROWS_NOTHING( manager.setProperty("AProp",5) )
    TS_ASSERT_THROWS( manager.setProperty("wefhui",5), Exception::NotFoundError )
    TS_ASSERT_THROWS( manager.setProperty("APROP",5.55), std::invalid_argument )
    TS_ASSERT_THROWS( manager.setProperty("APROP","value"), std::invalid_argument )
    TS_ASSERT_THROWS_NOTHING( manager.setProperty("AProp",1) )
  }

  void testExistsProperty()
  {
    Property *p = new PropertyWithValue<int>("sjfudh",0);
    TS_ASSERT( ! manager.existsProperty(p->name()) )
      Property *pp = new PropertyWithValue<double>("APROP",9.99);
    // Note that although the name of the property is the same, the type is different - yet it passes
    TS_ASSERT( manager.existsProperty(pp->name()) )
      delete p;
    delete pp;
  }

  void testValidateProperties()
  {
    TS_ASSERT( manager.validateProperties() )

      PropertyManagerHelper mgr;
    mgr.declareProperty("someProp","", new MandatoryValidator<std::string>);
    TS_ASSERT( ! mgr.validateProperties() )
  }

  void testGetPropertyValue()
  {
    TS_ASSERT( ! manager.getPropertyValue("APROP").compare("1") )
    TS_ASSERT_THROWS( manager.getPropertyValue("sdfshdu"), Exception::NotFoundError )
  }

  void testGetProperty()
  {
    Property *p = manager.getProperty("APROP");
    TS_ASSERT( p )
    TS_ASSERT( ! p->name().compare("aProp") )
    TS_ASSERT( ! p->value().compare("1") )
    TS_ASSERT( ! p->documentation().compare("") )
    TS_ASSERT( typeid(int) == *p->type_info() )

    TS_ASSERT_THROWS( p = manager.getProperty("werhui"), Exception::NotFoundError )

    int i;
    TS_ASSERT_THROWS_NOTHING( i = manager.getProperty("aprop") )
    TS_ASSERT_EQUALS( i, 1 );
    TS_ASSERT_THROWS( double dd = manager.getProperty("aprop"), std::runtime_error )
    std::string s = manager.getProperty("aprop");
    TS_ASSERT( ! s.compare("1") )
    double d;
    TS_ASSERT_THROWS_NOTHING( d = manager.getProperty("anotherProp") )
    TS_ASSERT_EQUALS( d, 1.11 );
    TS_ASSERT_THROWS( int ii = manager.getProperty("anotherprop"), std::runtime_error )
    std::string ss = manager.getProperty("anotherprop");
    // Note that some versions of boost::lexical_cast > 1.34 give a string such as
    // 9.9900000000000002 rather than 9.99. Converting back to a double however does
    // still give the correct 9.99.

    TS_ASSERT_EQUALS( ss.substr(0,4), std::string("1.11") )

    // This works, but CANNOT at present declare the string on a separate line and then assign
    //               (as I did for the int & double above)
    std::string sss = manager.getProperty("yetanotherprop");
    TS_ASSERT( ! sss.compare("itsValue") )
  }

  void testGetProperties()
  {
    std::vector<Property*> props = manager.getProperties();
    TS_ASSERT( props.size() == 3 )
    Property *p = props[0];
    TS_ASSERT( ! p->name().compare("aProp") )
    TS_ASSERT( ! p->value().compare("1") )
  }

private:
  PropertyManagerHelper manager;

};

#endif /*PROPERTYMANAGERTEST_H_*/
