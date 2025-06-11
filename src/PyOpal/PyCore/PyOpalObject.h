//
// Python API for PolynomialCoefficient (part of the multidimensional polynomial fitting routines)
//
// Copyright (c) 2008-2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef PyOpalObject_H
#define PyOpalObject_H

#include <Python.h>
#include <structmember.h>

#include <memory>
#include <exception>
#include <iostream>
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include <boost/mpl/front.hpp>

#include "Utilities/OpalException.h"
#include "Elements/OpalElement.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbsBeamline/Component.h"
#include "Attributes/Attributes.h"
#include "OpalConfigure/Configure.h"

extern Inform *gmsg;

namespace PyOpal {

/** PyOpalObjectNS namespace contains PyOpalObject, a wrapper for Object objects,
 *  and various supporting objects.
 *
 *  PyOpalObject<C>: basic element wrapper for C, which should be a subclass of
 *                OpalElement.
 *  AttributeType: enumeration of Opal Attribute Types (real, string, etc)
 *  AttributeDef: struct containing all of the things PyOpalObject needs to know
 *                about each attribute that should be exposed to the python api.
 *  PyElementGetProperty: call policy to handle access of an Attribute for a 
 *                python property
 *  PyElementSetProperty: call policy to handle setting of an Attribute for a
 *                python property
 *
 *  To wrap an OpalElement, say MyOpalObjectClass, you need to:
 *  1. Define static member data for the PyOpalObject<MyOpalObjectClass>. This 
 *     includes defining the attributes and setting a few options for different
 *     methods to expose.
 *  2. In BOOST_PYTHON_MODULE(my_opal_element_module), call 
 *
 *          PyOpalObject<MyOpalObjectType>.make_class()
 *
 *  Nb: apologies, this is heavy template stuff so almost everything has to go
 *  in the header file.
 */
namespace PyOpalObjectNS {

// forward declarations
template <class C> class PyOpalObject;
template <class C> struct PyOpalObjectGetProperty;
template <class C> struct PyOpalObjectSetProperty;

/** AttributeType is used to control conversion from python to OpalAttribute 
 *  - Float will convert to RealAttribute
 *  - String will convert to StringAttribute
 *  - Bool will convert to BoolAttribute (tho in python Bool is alias to long)
 *  - Long will convert to RealAttribute
 *  - VectorDouble will convert list to RealArray
 */
enum AttributeType {STRING, PREDEFINED_STRING, UPPER_CASE_STRING, STRING_LIST,
                    DOUBLE, BOOL, INT, FLOAT_LIST};

/** Maps the AttributeType to a string representation for docstrings/etc */
extern std::map<AttributeType, std::string> attributeName; // defined in PyOpalObject.cpp

/** AttributeDef defines an attribute
 *  opalName_m: the name of the opal Attribute
 *  pyName_m: the name that will be visible to the user in python. Properties
 *        should be_lower_case_with_underscores to comply with python API
 *  docString_m: docstring. If left empty (""), PyElement will generate a
 *        docstring like "py_name (type): Attribute help string"
 *  type_m: python type.
 */
struct AttributeDef {
    std::string opalName_m;
    std::string pyName_m;
    std::string docString_m;
    AttributeType type_m;
    // could also add a "read only" flag as some attributes are read only
};

/** PyOpalObject<C> is the basic wrapper class for Opal Objects
 *
 *  PyOpalObject<C> is a wrapper for an Opal Object C. Opal Attributes are 
 *  implemented in python as properties (public member data). Some default
 *  method calls can be enabled in the concrete implementation by setting flags.
 *  Docstring can be enabled by setting classDocString.
 *
 *  element_m: pointer to the base element
 *  attributes: vector of attributes. The order of the vector is the order that
 *      attributes will appear in the docstring
 *  classDocString: overall docstring. docstrings on each property and method 
 *      are added automatically by the framework.
  *
 *  property access is set up in make_class using add_property from
 *  boost::python. Normally, the "get" method can't take any data; this would
 *  mean that we can't pass an Attribute name for the OpalElement lookup. In
 *  order to get around this limitation, we have a "dummy" getter dummyGet.
 *  Then there is a "real" getter, which is hidden in the call policy for the
 *  property (a sort of pre- and post- decorator). We use PyOpalObjectGetProperty
 *  to hold the Attribute name and then dummyGet just tells PyElementGetProperty
 *  the address of the PyOpalObject. PyElementGetProperty does a postcall action to
 *  overwrite the return value with the "real" return value, based on the stored
 *  Attribute name. 
 *  
 *  Setters work in exactly the same way, except now we have dummySet and
 *  PyElementSetProperty to handle the interface.
 *  
 *  The routine then looks like:
 *     python::boost calls dummyGet;
 *           dummyGet tells PyElementGetProperty pointer to PyOpalObject
 *     python::boost calls PyElementGetProperty.postcall
 *           PyElementGetProperty.postcall calls PyOpalObject->getAttribute; where
 *           getAttribute takes the name of the Attribute as a string.
 *  It's a bit of a faff.
 */
template <class C>
class PyOpalObject {
public:
    typedef PyOpalObject<C> PyC; // for convenience
    /** Default constructor */
    inline PyOpalObject();
    /** Constructor taking the element as argument */
    PyOpalObject(std::shared_ptr<C> object) : object_m(object) {}
    /** Copying is disabled */
    PyOpalObject(const PyOpalObject<C>& rhs);
    /** Default destructor */
    virtual ~PyOpalObject() = default;

    /** This is the basic method to make a class for OpalObjects. It should
     *  normally be called when the module is declared, BOOST_PYTHON_MODULE.
     *
     *  Note that while in principle this could be static, in the end we need an
     *  instance of C so that we can access the docString from the attribute.
     */
    inline boost::python::class_<PyC> make_class(const char* className);

    /** This is the basic method to make a class *for elements*.
     * 
     *  This differs from make_class(...) because the setOpalName method sets
     *  the underlying element name at the same time.
     */
    inline boost::python::class_<PyC> make_element_class(const char* className);


    /** This is the basic method to make a class for generic Opal object type.
     * 
     *  In this one, the setOpalName method is not assigned, and user has to
     *  assign it (if they want it)
     */
    inline boost::python::class_<PyC> make_generic_class(const char* className);


    /** Add attributes to the python class.
     */
    template <class PYCLASS>
    void addAttributes(PYCLASS& pyclass);

    /** Add an "execute" method to the python class (e.g. if it is an ACTION)
     */
    template <class PYCLASS>
    void addExecute(PYCLASS& pyclass);

    /** Add a "register" method to the python class (to register against opal
     *  global data object)
     */
    template <class PYCLASS>
    void addRegister(PYCLASS& pyclass);

    /** Add a get_opal_name method to the python class (to get the opal internal
     *  string that uniquely identifies the object)
     */
    template <class PYCLASS>
    void addGetOpalName(PYCLASS& pyclass);

    /** Add a set_attributes method to the python class (to setup opalobject
     *  from a dictionary of kwargs)
     */
    template <class PYCLASS>
    void addSetAttributes(PYCLASS& pyclass);

    /** Add a set_opal_name method to the python class (to set the opal internal
     *  string that uniquely identifies the object)
     */
    template <class PYCLASS>
    void addSetOpalName(PYCLASS& pyclass);

    /** Add a set_opal_name method to the python class (to set the opal internal
     *  string that uniquely identifies the object). This version has bindings
     *  specifically for OpalElements
     */
    template <class PYCLASS>
    void addSetOpalElementName(PYCLASS& pyclass);

    /** Add a "get_opal_element" method to the python class (to overload as a
     *  PyOpalElement)
     * 
     *  This is required for OpalElement objects, in order that they can be used
     *  in a OPAL line object (for setting up beamlines having multiple
     *  elements). In general, if C is an OpalElement, developer should call
     *  addGetOpalElement in the module setup function.
     */
    template <class PYCLASS>
    void addGetOpalElement(PYCLASS& pyclass);

    /** Add a "get_field_value" method to the python class (for elements that
     *  expose a field)
     *
     *  - PYCLASS: instance of the class to which the field is added
     *  - distanceUnits: set the distance units so that UI is in OPAL units
     *  - timeUnits: set the time units so that UI is in OPAL units
     *  - bfieldUnits: set the bfield units so that UI is in OPAL units
     *  - efieldUnits: set the efield units so that UI is in OPAL units
     * 
     *  It is encouraged to expose the field map for field-type objects, for
     *  user to check that the field is as expected.
     *
     *  Note that conversion factors for units in distance, time, bfield and
     *  efield should be provided - the OPAL unit system should be exposed to
     *  user (see OPAL docs for current recommendation). Opal will call
     *  element->apply(R*distanceUnits, P, t*timeUnits, Efield, Bfield) and
     *  return Efield*efieldUnits, BField*bfield units.
     */
    template <class PYCLASS>
    void addGetFieldValue(PYCLASS& pyclass,
                          double distanceUnits, double timeUnits,
                          double bfieldUnits, double efieldUnits);

    /** dummyGet sets the object ptr for PyOpalObjectGetProperty but doesn't
     *  actually do the get(...)
     */
    template <class ValueType>
    ValueType dummyGet() const {PyOpalObjectGetProperty<C>::setObject(this); return ValueType();}

    /** Overload the method to perform any additional setup that can be made 
     *  before the OpalObject is accessed, after OpalObject->update() is called.
     */
    virtual void doSetup() {}

    /** dummySet sets the element ptr for PyOpalObjectSetProperty but doesn't
     *  actually do the set(...)
     */
    template <class ValueType>
    inline void dummySet(ValueType test);

    /** Get the value of an attribute
     *    - type: the type expected for python
     *    - opalName: the name of the Opal Attribute
     *  Returns a PyObject holding data; INCREF is called on the object i.e. it
     *  is a new object, not borrowed.
     *
     *  It is an error if the Attribute does not exist (throws OpalException)
     */
    PyObject* getAttribute(AttributeType type, std::string opalName) const;

    /** Set the value of an attribute
     *    - type: the type expected for python
     *    - opalName: the name of the Opal Attribute
     *    - value: Python object holding the value to be stored. value refcnt is
     *      not changed - caller is responsible for managing the refcnt.
     *  It is an error if the Attribute does not exist (throws OpalException)
     */
    void setAttribute(AttributeType type, std::string opalName, PyObject* value);

    /** Returns the Opal Object from the PyOpalObject */
    std::shared_ptr<C> getOpalShared() {return object_m;}

    /** Returns the Opal Object from the PyOpalObject (const version) */
    std::shared_ptr<C> getOpalShared() const {return object_m;}


protected:
    static std::vector<AttributeDef> attributes; /** class data (attributes) */
    static std::map<std::string, AttributeDef> pyNameToAttribute; /** maps python name to AttributeDef for fast lookup */
    static std::string classDocstring; /** class docstring */
    static bool converterRegistered; /** set to true if the converter has been registered */
    std::shared_ptr<C> object_m; /** pointer to the element */

    /** Generates a docstring from the attribute */
    std::string getDocString(AttributeDef& def);
    /** method to get field at a point in space-time */
    static boost::python::object getFieldValue(
                    PyOpalObjectNS::PyOpalObject<C>& pyobject,
                    double x, double y, double z, double t);
    /** method to get the opal name of an object */
    static std::string getOpalName(const PyOpalObject<C>& pyobject);
    /** method to set the opal name of an object */
    static void setOpalName(PyOpalObject<C>& pyobject, std::string name);
    /** method to set the opal name of an element - elements have two "names",
     *  one belonging to the OpalObject and one belonging to the Component */
    static void setOpalElementName(PyOpalObject<C>& pyobject, std::string name);
    /** method to call "execute" on an object */
    static void execute(PyOpalObject<C>& pyobject);
    /** method to call "register" on an object */
    static void registerObject(PyOpalObject<C>& pyobject);
    /** return the underlying element (used e.g. by PyLine) */
    static boost::python::object getPyOpalElement(PyOpalObject<C>& pyobject);
    /** set many attributes at once as a dict, stored in kwarg. 
     *  args should be exactly equal to (self,) */
    static boost::python::object setAttributes(boost::python::tuple args,
                                               boost::python::dict kwargs);
    // unit definitions for getFieldValue method
    static double distanceUnits_m;
    static double timeUnits_m;
    static double bfieldUnits_m;
    static double efieldUnits_m;
    static const std::string getFieldValueDocString;
};

/** apparently this is okay because it is a template and doesnt break one definition rule */
template <class C>
std::map<std::string, AttributeDef> PyOpalObject<C>::pyNameToAttribute;

/** Call update on a pyelement
 * 
 *  Note that C must be a subtype of OpalElement or some other class that has an
 *  C::update() method (i.e. not OpalObject).
 */
template <class C>
void update(PyOpalObjectNS::PyOpalObject<C> pyelement) {
    std::shared_ptr<C> objectPtr = pyelement.getOpalShared();
    objectPtr->update();
}

template <class C>
template <class ValueType>
void PyOpalObject<C>::dummySet(ValueType test) {
    (void)test; // disable Wunused-parameter
    PyOpalObjectSetProperty<C>::setObject(this);
}

template <class C>
boost::python::object PyOpalObject<C>::getFieldValue(
            PyOpalObjectNS::PyOpalObject<C>& pyobject, 
            double x, double y, double z, double t) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    pyobject.doSetup();
    ElementBase* element = objectPtr->getElement();
    Component* component = dynamic_cast<Component*>(element);
    if (component == nullptr) {
        throw OpalException("PyElement<C>::getFieldValue",
                            "Failed to deduce Component from ElementBase.");
    }
    Vector_t R(x*distanceUnits_m, y*distanceUnits_m, z*distanceUnits_m);
    Vector_t P(0.0, 0.0, 0.0);
    Vector_t B;
    Vector_t E;
    t *= timeUnits_m;
    bool outOfBounds = component->apply(R, P, t, E, B);
    return boost::python::make_tuple(outOfBounds,
                    B[0]*bfieldUnits_m, B[1]*bfieldUnits_m, B[2]*bfieldUnits_m,
                    E[0]*efieldUnits_m, E[1]*efieldUnits_m, E[2]*efieldUnits_m);
}

/** Pass a dictionary of attributes to the pyobject using the python *args **kwargs semantics */
template <class C>
boost::python::object PyOpalObject<C>::setAttributes(boost::python::tuple args,
                                                     boost::python::dict kwargs) {
    if (boost::python::len(args) == 0) {
        throw OpalException("PyOpalObject::setAttributes",
            "Did not find class instance (i.e. 'self') when calling "
            "set_attributes");
    } else if (boost::python::len(args) > 1) {
        throw OpalException("PyOpalObject::setAttributes",
            "set_attributes cannot take any non-keyword args (except the class"
            " instance");
    }
    PyOpalObject<C> self = boost::python::extract<PyOpalObject<C> >(args[0]);
    boost::python::list key_list = kwargs.keys();
    for (boost::python::ssize_t i = 0; i < boost::python::len(key_list); ++i) {
        boost::python::object key = key_list[i];
        std::string keyStr = boost::python::extract<std::string>(key);
        boost::python::object value = kwargs[key];
        if (PyOpalObject<C>::pyNameToAttribute.find(keyStr) ==
                                    PyOpalObject<C>::pyNameToAttribute.end()) {
            throw OpalException("PyOpalObject::setAttributes",
                                "Did not recognise attribute '"+keyStr+"'");
        }
        AttributeDef att = PyOpalObject<C>::pyNameToAttribute[keyStr];
        self.setAttribute(att.type_m, att.opalName_m, value.ptr());
    }
    boost::python::object obj;
    return obj;
}

/** Helper class to handle getting Attributes from python */
template <class C>
struct PyOpalObjectGetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyOpalObjectGetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyOpalObjectGetProperty() {;}

    /** postcall action
     *  - ArgumentPackage: this is really a PyObject* holding the argument calls
     *  - result: filled by the postcall operation.
     *  Sets the element_m to nullptr. Note that the whole element_m stuff is a bit
     *  hacky. Probably possible to extract the element from ArgumentPackage, I
     *  didn't manage to do it yet.
     */
    template <class ArgumentPackage>
    PyObject* postcall(ArgumentPackage const&, PyObject* result);

    /** Set pointer to the element; should be called before each postcall*/
    static void setObject(const PyOpalObject<C>* object) {object_m = object;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static const PyOpalObject<C>* object_m;
};

/** Helper class to handle setting Attributes from python */
template <class C>
struct PyOpalObjectSetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyOpalObjectSetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyOpalObjectSetProperty() {;}

    /** postcall action
     *  - ArgumentPackage: this is really a PyObject* holding the argument calls
     *  - result: filled by the postcall operation (with PyNone).
     *  Sets the element_m to nullptr. Note that the whole element_m stuff is a bit
     *  hacky. Probably possible to extract the element from ArgumentPackage, I
     *  didn't manage to do it yet.
     */
    template <class ArgumentPackage>
    PyObject* postcall(ArgumentPackage const& args, PyObject* result);

    /** Set pointer to the element; should be called before each postcall*/
    static void setObject(PyOpalObject<C>* object) {object_m = object;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static PyOpalObject<C>* object_m;
};

///////// templated implementations for PyElement ////////////////


template <class C>
void PyOpalObject<C>::execute(PyOpalObject<C>& pyobject) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->execute();
}

template <class C>
void PyOpalObject<C>::registerObject(PyOpalObjectNS::PyOpalObject<C>& pyobject) {
    C* wrappedC = pyobject.getOpalShared().get();
    Object* objectPtr = dynamic_cast<Object*>(wrappedC);
    if (objectPtr == nullptr) {
        throw OpalException("PyOpalObject<C>::registerObject",
                            "Trying to register something that was not a Opal Object");
    }
    //Object* objectPtr = &(*pyobject.getOpalShared());
    OpalData::getInstance()->define(objectPtr);
}

template <class C>
std::string PyOpalObject<C>::getOpalName(const PyOpalObject<C>& pyobject) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    return objectPtr->getOpalName();
}

template <class C>
void PyOpalObject<C>::setOpalName(PyOpalObject<C>& pyobject, std::string name) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->setOpalName(name);
}

template <class C>
void PyOpalObject<C>::setOpalElementName(PyOpalObject<C>& pyobject, std::string name) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->setOpalName(name);
    if (objectPtr->getElement()) { // some OpalThings don't define an Element
        objectPtr->getElement()->setName(name);
    }
}

template <class C>
boost::python::object PyOpalObject<C>::getPyOpalElement(PyOpalObjectNS::PyOpalObject<C>& pyobject) {
    std::shared_ptr<OpalElement> elementPtr =
            std::dynamic_pointer_cast<OpalElement, C>(pyobject.getOpalShared());
    if (elementPtr.get() == nullptr) {
        throw OpalException("PyOpalObject<C>::getPyOpalElement",
                            "Wrapped object was not an OpalElement");
    }
    PyOpalObject<OpalElement> element(elementPtr);
    boost::python::object pyelement(element);
    return pyelement;
}


// defined in PyOpalElement
template <>
PyOpalObject<OpalElement>::PyOpalObject();

template <class C>
PyOpalObject<C>::PyOpalObject() : object_m(new C) {}

template <class C>
PyObject* PyOpalObject<C>::getAttribute(AttributeType type, std::string opalName) const {
    if (!object_m) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
                            "Object was not initialised");       
    }
    Attribute* attribute = object_m->findAttribute(opalName);
    if (attribute == nullptr) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
                            "Failed to parse attribute "+opalName);
    }
    PyObject* pyvalue;
    // I spent quite a bit of time trying to get this type unwinding to work
    // using templates (so it is done at compile time). In the end I couldn't
    // fight the template syntax and had to do it at runtime using an enum; it's
    // not so bad - the memory footprint is smaller and, tbh, template syntax is
    // horrible so might be easier to use.
    switch (type) {
        case DOUBLE: 
        {
            double value = Attributes::getReal(*attribute);
            pyvalue = PyFloat_FromDouble(value);
            break;
        }
        case INT:
        {
            double value = Attributes::getReal(*attribute);
            pyvalue = PyLong_FromDouble(value);
            break;
        }
        case STRING:
        case PREDEFINED_STRING:
        case UPPER_CASE_STRING:
        {
            std::string value = Attributes::getString(*attribute);
            pyvalue = PyUnicode_FromString(value.c_str());
            break;
        }
        case BOOL:
        {
            bool value = Attributes::getBool(*attribute);
            if (value) {
                pyvalue = Py_True;
            } else {
                pyvalue = Py_False;        
            }
            break;
        }
        case FLOAT_LIST:
        {
            std::vector<double> value = Attributes::getRealArray(*attribute);
            pyvalue = PyList_New(value.size());
            for (size_t i = 0; i < value.size(); ++i) {
                PyList_SetItem(pyvalue, i, PyFloat_FromDouble(value[i]));
            }
            break;
        }
        case STRING_LIST:
        {
            std::vector<std::string> value = Attributes::getStringArray(*attribute);
            pyvalue = PyList_New(value.size());
            for (size_t i = 0; i < value.size(); ++i) {
                PyList_SetItem(pyvalue, i, PyUnicode_FromString(value[i].c_str()));
            }
            break;
        }
        default:
        {
            throw OpalException("PyOpalObject<C>::getAttribute",
                            "Attribute type "+attributeName[type]+" not implemented");
        }
    }
    Py_INCREF(pyvalue);
    return pyvalue;
}


template <class C>
void PyOpalObject<C>::setAttribute(AttributeType type, std::string opalName, PyObject* pyvalue) {
    if (!object_m) {
        throw OpalException("PyOpalObject<C>::setAttribute",
                            "Element was not initialised");
    }
    Attribute* attribute = object_m->findAttribute(opalName);
    if (attribute == nullptr) {
        throw OpalException("PyOpalObject<C>::setAttribute",
                            "Failed to parse attribute "+opalName);
    }
    switch (type) {
        case DOUBLE:
        {
            double value = PyFloat_AsDouble(pyvalue);
            Attributes::setReal(*attribute, value);
            break;
        }
        case INT:
        {
            double value = PyLong_AsDouble(pyvalue);
            Attributes::setReal(*attribute, value);
            break;
        }
        case STRING:
        {
            std::string value = PyUnicode_AsUTF8(pyvalue);
            Attributes::setString(*attribute, value);
            break;
        }
        case PREDEFINED_STRING:
        {
            // predefined string is a string with a list of options
            // note that value checking never happens - method exists in 
            // src/Attributes/PredefinedString.h but requires a Statement object
            // which seems hard to construct (Classic/Parser/Statement.h)
            std::string value = PyUnicode_AsUTF8(pyvalue);
            Attributes::setPredefinedString(*attribute, value);
            break;
        }
        case UPPER_CASE_STRING:
        {
            std::string value = PyUnicode_AsUTF8(pyvalue);
            Attributes::setUpperCaseString(*attribute, value);
            break;
        }
        case BOOL:
        {
            bool value = PyObject_IsTrue(pyvalue);
            Attributes::setBool(*attribute, value);
            break;
        }
        case FLOAT_LIST:
        {
            Py_ssize_t listSize = PyList_Size(pyvalue);
            std::vector<double> value(listSize);
            for (Py_ssize_t i = 0; i < listSize; ++i) {
                double value_i = PyFloat_AsDouble(PyList_GetItem(pyvalue, i));
                value[i] = value_i;
            }
            Attributes::setRealArray(*attribute, value);
            break;
        }
        case STRING_LIST:
        {
            Py_ssize_t listSize = PyList_Size(pyvalue);
            std::vector<std::string> value(listSize);
            for (Py_ssize_t i = 0; i < listSize; ++i) {
                PyObject* pyvalue_i = PyList_GetItem(pyvalue, i);
                std::string value_i = PyUnicode_AsUTF8(pyvalue_i);
                value[i] = value_i;
            }
            Attributes::setStringArray(*attribute, value);
            break;
        }
        default:
        {
            throw OpalException("PyOpalObject<C>::setAttribute",
                                "Attribute type "+attributeName[type]+" not implemented");
        }
    }
}

template <class C>
PyOpalObject<C>::PyOpalObject(const PyOpalObject<C>& rhs) : object_m(rhs.object_m) { 
}

template <class C>
std::string PyOpalObject<C>::getDocString(AttributeDef& def) {
    Attribute* attribute = object_m->findAttribute(def.opalName_m);
    if (attribute == nullptr) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
                            "Failed to parse attribute "+def.opalName_m);
    }
    std::string docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+attribute->getHelp();
    if (def.docString_m != "") {
        docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+def.docString_m;
    }
    return docString;
}

template <class C>
boost::python::class_<PyOpalObject<C> > PyOpalObject<C>::make_generic_class(const char* className) {
    // WARNING - boost::python is bugged so that in module initialisation, 
    // errors are not handled correctly
    //    https://github.com/boostorg/python/issues/280
    // so if you get the attribute name wrong you get a bad error
    // hence I do some local error handling
    typedef boost::python::class_<PyOpalObject<C> > PyClass;
    boost::python::docstring_options docop(true, true, false); // user_def, py_sig, cpp_sig
    PyClass pyclass = PyClass(className);
    try {
        addAttributes(pyclass);
        addGetOpalName(pyclass);
        addSetAttributes(pyclass);
    } catch (OpalException& exc) {
        std::cerr << "Failed to initialise class because '" << exc.what() 
                  << "'" << std::endl;
        throw exc;
    } catch (std::exception& exc) {
        std::cerr << "Failed to initialise class because '" << exc.what() 
                  << "'" << std::endl;
        throw exc;
    }
    if (gmsg == nullptr) {
        std::cerr <<
            "Failed to initialise class because gmsg was not initialised.\n"
            "PyOpal::Globals::Initialise() must be called before attempting to\n"
            "initialise classes" << std::endl;
        throw OpalException("PyOpalObject::make_generic_class",
                            "Globals not initialised");
    }
    return pyclass;
}

template <class C>
boost::python::class_<PyOpalObject<C> > PyOpalObject<C>::make_class(const char* className) {
    auto pyclass = make_generic_class(className);
    addSetOpalName(pyclass);
    return pyclass;
}

template <class C>
boost::python::class_<PyOpalObject<C> > PyOpalObject<C>::make_element_class(const char* className) {
    auto pyclass = make_generic_class(className);
    addSetOpalElementName(pyclass);
    addGetOpalElement(pyclass);
    return pyclass;
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addExecute(PYCLASS& pyclass) {
    pyclass.def("execute", &PyOpalObject<C>::execute);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addRegister(PYCLASS& pyclass) {
    pyclass.def("register", &PyOpalObject<C>::registerObject);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addGetOpalName(PYCLASS& pyclass) {
    pyclass.def("get_opal_name", &PyOpalObject<C>::getOpalName);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addSetAttributes(PYCLASS& pyclass) {
    pyclass.def("set_attributes",
                boost::python::raw_function(&PyOpalObject<C>::setAttributes, 0)
                );
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addSetOpalName(PYCLASS& pyclass) {
    pyclass.def("set_opal_name", &PyOpalObject<C>::setOpalName);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addSetOpalElementName(PYCLASS& pyclass) {
    pyclass.def("set_opal_name", &PyOpalObject<C>::setOpalElementName);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addGetOpalElement(PYCLASS& pyclass) {
    pyclass.def("get_opal_element", &PyOpalObject<C>::getPyOpalElement);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addGetFieldValue(PYCLASS& pyclass,
                                      double distanceUnits, double timeUnits,
                                      double efieldUnits, double bfieldUnits) {
    distanceUnits_m = distanceUnits;
    timeUnits_m = timeUnits;
    bfieldUnits_m = bfieldUnits;
    efieldUnits_m = efieldUnits;
    pyclass.def("get_field_value",
                &PyOpalObject<C>::getFieldValue,
                boost::python::args("x", "y", "z", "t"),
                getFieldValueDocString.c_str());
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addAttributes(PYCLASS& pyclass) {
    for (std::vector<AttributeDef>::iterator iter = attributes.begin(); iter != attributes.end(); ++iter) {
        pyNameToAttribute[iter->pyName_m] = *iter;
        PyOpalObjectGetProperty<C> getProp(iter->type_m, iter->opalName_m);
        PyOpalObjectSetProperty<C> setProp(iter->type_m, iter->opalName_m);
        std::string docString = getDocString(*iter);
        std::string pyname = iter->pyName_m.c_str();
        switch (iter->type_m) {
            case DOUBLE:
            {
                pyclass.add_property(pyname.c_str(),
                                     boost::python::make_function(&PyC::dummyGet<double>, getProp),
                                     boost::python::make_function(&PyC::dummySet<double>, setProp),
                                     docString.c_str()
                );
                break;
            }
            case INT:
            {
                pyclass.add_property(pyname.c_str(),
                                     boost::python::make_function(&PyC::dummyGet<int>, getProp),
                                     boost::python::make_function(&PyC::dummySet<int>, setProp),
                                     docString.c_str()
                );
                break;
            }
            case STRING:
            case PREDEFINED_STRING:
            case UPPER_CASE_STRING:
            {
                pyclass.add_property(pyname.c_str(),
                                     boost::python::make_function(&PyC::dummyGet<std::string>, getProp),
                                     boost::python::make_function(&PyC::dummySet<std::string>, setProp),
                                     docString.c_str()
                );
                break;
            }
            case BOOL:
            {
                pyclass.add_property(pyname.c_str(),
                                     boost::python::make_function(&PyC::dummyGet<bool>, getProp),
                                     boost::python::make_function(&PyC::dummySet<bool>, setProp),
                                     docString.c_str()
                );
                break;
            }
            case FLOAT_LIST:
            case STRING_LIST:
            {
                pyclass.add_property(pyname.c_str(),
                                     boost::python::make_function(&PyC::dummyGet<boost::python::list>, getProp),
                                     boost::python::make_function(&PyC::dummySet<boost::python::list>, setProp),
                                     docString.c_str()
                );
                break;
            } 
            default:
            {
                // Looks like exception handling doesn't work properly at module
                // import time so this may not be handled politely - thrown as an 
                // unrecognised SystemError
                throw OpalException("PyOpalObject<C>::addAttributes", "Type not implemented");
            }
        }
    }
}

template <class C> double PyOpalObject<C>::distanceUnits_m = 1;
template <class C> double PyOpalObject<C>::timeUnits_m = 1;
template <class C> double PyOpalObject<C>::bfieldUnits_m = 1;
template <class C> double PyOpalObject<C>::efieldUnits_m = 1;

template <class C>  const std::string PyOpalObject<C>::getFieldValueDocString =
  "Get the field value at a point in the field map.\n"
  "\n"
  "The field lookup is performed against the last RINGDEFINITION that was\n"
  "instantiated. This should be instantiated by calling\n"
  "pyopal.parser.initialise_from_opal_file\n"
  "\n"
  "Parameters\n"
  "----------\n"
  "x : float\n"
  "    x position [m]\n"
  "y : float\n"
  "    y position [m]\n"
  "z : float\n"
  "    z position [m]\n"
  "t: float\n"
  "    time [ns]\n"
  "\n"
  "Returns\n"
  "-------\n"
  "The function returns a tuple containing 7 values:\n"
  "out of bounds : int\n"
  "    1 if the event was out of the field map boundary, else 0.\n"
  "Bx : float\n"
  "    x magnetic field [T]\n"
  "By : float\n"
  "    y magnetic field [T]\n"
  "Bz : float\n"
  "    z magnetic field [T]\n"
  "Ex : float\n"
  "    x electric field [MV/m]\n"
  "Ey : float\n"
  "    y electric field [MV/m]\n"
  "Ez : float\n"
  "    z electric field [MV/m]\n";


///////// templated implementations for PyOpalObjectGetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyOpalObjectGetProperty<C>::postcall(ArgumentPackage const&, PyObject* result) {
        Py_DECREF(result);
        result = object_m->getAttribute(type_m, opalName_m);
        return result;
}

template <class C>
const PyOpalObject<C>* PyOpalObjectGetProperty<C>::object_m = nullptr;

///////// templated implementations for PyOpalObjectSetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyOpalObjectSetProperty<C>::postcall(ArgumentPackage const& args, PyObject* result) {
    PyObject* value;
    PyObject* pyObject; // this is a direct pointer to the C but I don't know how to unwrap it...
    if (!PyArg_ParseTuple(args, "OO", &pyObject, &value)) {
        return nullptr; // ParseTuple sets the error message
    }
    Py_DECREF(result);
    object_m->setAttribute(type_m, opalName_m, value);
    object_m = nullptr;
    Py_RETURN_NONE;
}

template <class C>
PyOpalObject<C>* PyOpalObjectSetProperty<C>::object_m = nullptr;

} // PyOpalObject
} // PyOpal

#endif // PyOpalObject_H
