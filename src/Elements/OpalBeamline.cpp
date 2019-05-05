#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Bend.h"
#include "Structure/MeshGenerator.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <fstream>
using namespace std;

OpalBeamline::OpalBeamline():
    elements_m(),
    prepared_m(false),
    containsSource_m(false)
{
}

OpalBeamline::OpalBeamline(const Vector_t& origin,
                           const Quaternion& rotation):
    elements_m(),
    prepared_m(false),
    containsSource_m(false),
    coordTransformationTo_m(origin, rotation)
{
}

OpalBeamline::~OpalBeamline() {
    elements_m.clear();
}

CompVec OpalBeamline::dummy_list_m = CompVec();
OpalSection OpalBeamline::dummy_section_m = OpalSection(dummy_list_m, 0., 0.);

std::set<std::shared_ptr<Component>> OpalBeamline::getElements(const Vector_t &x) {
    std::set<std::shared_ptr<Component> > elementSet;
    FieldList::iterator it = elements_m.begin();
    const FieldList::iterator end = elements_m.end();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        Vector_t r = (*it).getCoordTransformationTo().transformTo(x);

        if (element->isInside(r)) {
            elementSet.insert(element);
        }
    }

    return elementSet;
}

void OpalBeamline::getKFactors(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &KR, Vector_t &KT) {

}

unsigned long OpalBeamline::getFieldAt(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &E, Vector_t &B) {

    unsigned long rtv = 0x00;

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(const Vector_t &position,
                                       const Vector_t &momentum,
                                       const double &t,
                                       Vector_t &Ef,
                                       Vector_t &Bf) {
    unsigned long rtv = 0x00;

    std::set<std::shared_ptr<Component>> elements = getElements(position);

    std::set<std::shared_ptr<Component>>::const_iterator it = elements.begin();
    const std::set<std::shared_ptr<Component>>::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        ElementBase::ElementType type = (*it)->getType();
        if (type == ElementBase::MONITOR ||
            type == ElementBase::MARKER ||
            type == ElementBase::CCOLLIMATOR ||
            type == ElementBase::DIAGNOSTIC) continue;

        Vector_t localR = transformToLocalCS(*it, position);
        Vector_t localP = rotateToLocalCS(*it, momentum);
        Vector_t localE(0.0), localB(0.0);

        (*it)->applyToReferenceParticle(localR, localP, t, localE, localB);

        Ef += rotateFromLocalCS(*it, localE);
        Bf += rotateFromLocalCS(*it, localB);
    }

    //         if(section.hasWake()) {
    //             rtv |= BEAMLINE_WAKE;
    //         }
    //         if(section.hasParticleMatterInteraction()) {
    //             rtv |= BEAMLINE_PARTICLEMATTERINTERACTION;
    //         }

    return rtv;
}

void OpalBeamline::switchElements(const double &min, const double &max, const double &kineticEnergy, const bool &nomonitors) {

    FieldList::iterator fprev;
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        // don't set online monitors if the centroid of the bunch is allready inside monitor
        // or if explicitly not desired (eg during auto phasing)
        if(flit->getElement()->getType() == ElementBase::MONITOR) {
            double spos = (max + min) / 2.;
            if(!nomonitors && spos < (*flit).getStart()) {
                if(!(*flit).isOn() && max > (*flit).getStart()) {
                    (*flit).setOn(kineticEnergy);
                }
            }

        } else {
            if(!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
                (*flit).setOn(kineticEnergy);
            }

            //check if multiple degraders follow one another with no other elements in between
            //if element is off and it is a degrader
            if (!(*flit).isOn() && flit->getElement()->getType() == ElementBase::DEGRADER) {
                //check if previous element: is on, is a degrader, ends where new element starts
                if ((*fprev).isOn() && fprev->getElement()->getType() == ElementBase::DEGRADER &&
                    ((*fprev).getEnd() + 0.01 > (*flit).getStart()) ) {
                    (*flit).setOn(kineticEnergy);
                }
            }
        }

        fprev = flit;
    }
}


void OpalBeamline::switchAllElements() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if(!(*flit).isOn()) {
            (*flit).setOn(-1.0);
        }
    }

}

void OpalBeamline::switchElementsOff(const double &min, ElementBase::ElementType eltype) {
    if(eltype == ElementBase::ANY) {
        for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if((*flit).isOn() && min >= (*flit).getEnd()) {
                (*flit).setOff();
            }

        }
    } else {
        for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if((*flit).isOn() && min >= (*flit).getEnd() && (*flit).getElement()->getType() == eltype) {
                (*flit).setOff();
            }
        }
    }
}

void OpalBeamline::switchElementsOff() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit)
        (*flit).setOff();
}

void OpalBeamline::prepareSections() {
    if (elements_m.size() == 0) {
        prepared_m = true;
        return;
    }
    prepared_m = true;
}

void OpalBeamline::print(Inform &msg) const {
}


double OpalBeamline::calcBeamlineLength() {
    return 0.0;
}

void OpalBeamline::swap(OpalBeamline & rhs) {
    std::swap(elements_m, rhs.elements_m);
    std::swap(prepared_m, rhs.prepared_m);
    std::swap(coordTransformationTo_m, rhs.coordTransformationTo_m);
}

void OpalBeamline::merge(OpalBeamline &rhs) {
    elements_m.insert(elements_m.end(),
                      rhs.elements_m.begin(),
                      rhs.elements_m.end());
    prepared_m = false;
    containsSource_m = containsSource_m || rhs.containsSource_m;
}


FieldList OpalBeamline::getElementByType(ElementBase::ElementType type) {
    if (type == ElementBase::ANY) {
        return elements_m;
    }

    FieldList elements_of_requested_type;
    for(FieldList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++ fit) {
        if((*fit).getElement()->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}

void OpalBeamline::compute3DLattice() {
    static unsigned int order = 0;
    FieldList::iterator it = elements_m.begin();
    const FieldList::iterator end = elements_m.end();

    unsigned int minOrder = order;
    {
        double endPriorPathLength = 0.0;
        CoordinateSystemTrafo currentCoordTrafo = coordTransformationTo_m;

        elements_m.sort([](const ClassicField& a, const ClassicField& b) {
                double edgeA = 0.0, edgeB = 0.0;
                if (a.getElement()->isElementPositionSet())
                    edgeA = a.getElement()->getElementPosition();

                if (b.getElement()->isElementPositionSet())
                    edgeB = b.getElement()->getElementPosition();

                return edgeA < edgeB;
            });
        it = elements_m.begin();
        for (; it != end; ++ it) {
            if ((*it).isPositioned()) {
                continue;
            }
            (*it).order_m = minOrder;
            std::shared_ptr<Component> element = (*it).getElement();
            if (element->getType() == ElementBase::SBEND ||
                element->getType() == ElementBase::RBEND ||
                element->getType() == ElementBase::RBEND3D) {

                double beginThisPathLength = element->getElementPosition();
                Vector_t beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);
                BendBase * bendElement = static_cast<BendBase*>(element.get());
                double thisLength = bendElement->getChordLength();
                double bendAngle = bendElement->getBendAngle();
                double entranceAngle = bendElement->getEntranceAngle();
                double arcLength = (thisLength * std::abs(bendAngle) / (2 * sin(std::abs(bendAngle) / 2)));

                double rotationAngleAboutZ = bendElement->getRotationAboutZ();
                Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                            sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

                Vector_t effectiveRotationAxis = rotationAboutZ.rotate(Vector_t(0, -1, 0));
                effectiveRotationAxis /= euclidean_norm(effectiveRotationAxis);

                Quaternion_t rotationAboutAxis(cos(0.5 * bendAngle),
                                               sin(0.5 * bendAngle) * effectiveRotationAxis);
                Quaternion_t halfRotationAboutAxis(cos(0.25 * bendAngle),
                                                   sin(0.25 * bendAngle) * effectiveRotationAxis);
                Quaternion_t entryFaceRotation(cos(0.5 * entranceAngle),
                                               sin(0.5 * entranceAngle) * effectiveRotationAxis);

                if (!Options::idealized) {
                    std::vector<Vector_t> truePath = bendElement->getDesignPath();
                    Quaternion_t directionExitHardEdge(cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                                                       sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                    Vector_t exitHardEdge = thisLength * directionExitHardEdge.rotate(Vector_t(0, 0, 1));
                    double distanceEntryHETruePath = euclidean_norm(truePath.front());
                    double distanceExitHETruePath = euclidean_norm(rotationAboutZ.rotate(truePath.back()) - exitHardEdge);
                    double pathLengthTruePath = (*it).getEnd() - (*it).getStart();
                    arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
                }

                Vector_t chord = thisLength * halfRotationAboutAxis.rotate(Vector_t(0, 0, 1));
                Vector_t endThis3D = beginThis3D + chord;
                double endThisPathLength = beginThisPathLength + arcLength;

                CoordinateSystemTrafo fromEndLastToBeginThis(beginThis3D,
                                                             (entryFaceRotation * rotationAboutZ).conjugate());
                CoordinateSystemTrafo fromEndLastToEndThis(endThis3D,
                                                           rotationAboutAxis.conjugate());

                (*it).setCoordTransformationTo(fromEndLastToBeginThis * currentCoordTrafo);

                currentCoordTrafo = (fromEndLastToEndThis * currentCoordTrafo);

                endPriorPathLength = endThisPathLength;
            }
        }
    }

    double endPriorPathLength = 0.0;
    CoordinateSystemTrafo currentCoordTrafo = coordTransformationTo_m;

    it = elements_m.begin();
    for (; it != end; ++ it) {
        if ((*it).isPositioned()) continue;

        (*it).order_m = order ++;

        std::shared_ptr<Component> element = (*it).getElement();
        double beginThisPathLength = element->getElementPosition();
        double thisLength = element->getElementLength();
        Vector_t beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);

        if (element->getType() == ElementBase::MONITOR) {
            beginThis3D(2) -= 0.5 * thisLength;
        }

        Vector_t endThis3D;
        if (element->getType() == ElementBase::SBEND ||
            element->getType() == ElementBase::RBEND ||
            element->getType() == ElementBase::RBEND3D) {

            BendBase * bendElement = static_cast<BendBase*>(element.get());
            thisLength = bendElement->getChordLength();
            double bendAngle = bendElement->getBendAngle();

            double rotationAngleAboutZ = bendElement->getRotationAboutZ();
            Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                        sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

            Vector_t effectiveRotationAxis = rotationAboutZ.rotate(Vector_t(0, -1, 0));
            effectiveRotationAxis /= euclidean_norm(effectiveRotationAxis);

            Quaternion_t rotationAboutAxis(cos(0.5 * bendAngle),
                                           sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion halfRotationAboutAxis(cos(0.25 * bendAngle),
                                             sin(0.25 * bendAngle) * effectiveRotationAxis);

            double arcLength = (thisLength * std::abs(bendAngle) /
                                (2 * sin(bendAngle / 2)));
            if (!Options::idealized) {
                std::vector<Vector_t> truePath = bendElement->getDesignPath();
                double entranceAngle = bendElement->getEntranceAngle();
                Quaternion_t directionExitHardEdge(cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                                                   sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                Vector_t exitHardEdge = thisLength * directionExitHardEdge.rotate(Vector_t(0, 0, 1));
                double distanceEntryHETruePath = euclidean_norm(truePath.front());
                double distanceExitHETruePath = euclidean_norm(rotationAboutZ.rotate(truePath.back()) - exitHardEdge);
                double pathLengthTruePath = (*it).getEnd() - (*it).getStart();
                arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
            }

            endThis3D = (beginThis3D +
                         halfRotationAboutAxis.rotate(Vector_t(0, 0, thisLength)));
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D,
                                                       rotationAboutAxis.conjugate());
            currentCoordTrafo = fromEndLastToEndThis * currentCoordTrafo;

            endPriorPathLength = beginThisPathLength + arcLength;
        } else {
            double rotationAngleAboutZ = (*it).getElement()->getRotationAboutZ();
            Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                        sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

            CoordinateSystemTrafo fromLastToThis(beginThis3D, rotationAboutZ);

            (*it).setCoordTransformationTo(fromLastToThis * currentCoordTrafo);
        }

        (*it).fixPosition();
    }

    elements_m.sort(ClassicField::SortAsc);
}

void OpalBeamline::plot3DLattice() {
    auto opal = OpalData::getInstance();
    if (opal->isOptimizerRun() || Ippl::myNode() != 0) return;

    elements_m.sort([](const ClassicField& a, const ClassicField& b) {
            double edgeA = 0.0, edgeB = 0.0;
            if (a.getElement()->isElementPositionSet())
                edgeA = a.getElement()->getElementPosition();

            if (b.getElement()->isElementPositionSet())
                edgeB = b.getElement()->getElementPosition();

            return edgeA < edgeB;
        });

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    Quaternion rotDiagonal(0.5, 0.5 * Vector_t(-1, 1, -1));

    Vector_t origin = rotDiagonal.rotate(coordTransformationTo_m.getOrigin());
    Vector_t direction = rotDiagonal.rotate(coordTransformationTo_m.rotateFrom(Vector_t(0, 0, 1)));
    Vector_t minX = Vector_t(999999.9), maxX(-999999.9);
    std::map<std::string, std::vector<Vector_t > > elementCorners;

    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        CoordinateSystemTrafo toBegin = (*it).getCoordTransformationTo();
        std::vector<Vector_t> corners;

        if (element->getType() == ElementBase::RBEND || element->getType() == ElementBase::SBEND) {
            std::vector<Vector_t> outline = static_cast<const Bend*>(element.get())->getOutline();

            for (auto point: outline) {
                corners.push_back(rotDiagonal.rotate(toBegin.transformFrom(point)));
            }
        } else {
            CoordinateSystemTrafo toEnd = element->getEdgeToEnd() * toBegin;
            auto aperture = element->getAperture();
            double elementHeightFront = aperture.second[0];
            double elementHeightBack = aperture.second[0] * aperture.second[2];

            corners.push_back(rotDiagonal.rotate(toEnd.transformFrom(Vector_t(0.0))));
            corners.push_back(rotDiagonal.rotate(toBegin.transformFrom(Vector_t(0))));
            corners.push_back(rotDiagonal.rotate(toBegin.transformFrom(-elementHeightFront * Vector_t(1, 0, 0))));
            corners.push_back(rotDiagonal.rotate(toEnd.transformFrom(-elementHeightBack * Vector_t(1, 0, 0))));
            corners.push_back(rotDiagonal.rotate(toEnd.transformFrom(Vector_t(0.0))));
            corners.push_back(rotDiagonal.rotate(toEnd.transformFrom(elementHeightBack * Vector_t(1, 0, 0))));
            corners.push_back(rotDiagonal.rotate(toBegin.transformFrom(elementHeightFront * Vector_t(1, 0, 0))));
            corners.push_back(rotDiagonal.rotate(toBegin.transformFrom(Vector_t(0))));
        }

        elementCorners.insert(std::make_pair(element->getName(), corners));
        const unsigned int numCorners = corners.size();
        for (unsigned int i = 0 ; i < numCorners; ++ i) {
            const Vector_t & X = corners[i];

            if (X(0) < minX(0)) minX(0) = X(0);
            else if (X(0) > maxX(0)) maxX(0) = X(0);

            if (X(1) < minX(1)) minX(1) = X(1);
            else if (X(1) > maxX(1)) maxX(1) = X(1);
        }
    }

    std::ofstream gpl;
    std::string fileName = "data/" + OpalData::getInstance()->getInputBasename() + "_ElementPositions.gpl";
    if (Options::openMode == Options::APPEND && boost::filesystem::exists(fileName)) {
        gpl.open(fileName, std::ios_base::app);
    } else {
        gpl.open(fileName);
    }
    gpl.precision(8);

    it = elements_m.begin();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();

        if (element->getType() != ElementBase::DRIFT) {
            const std::vector<Vector_t> &corners = elementCorners[element->getName()];
            const unsigned int numCorners = corners.size();

            gpl << "# " << element->getName() << "\n";
            for (unsigned int i = 0; i < numCorners; ++ i) {
                gpl << std::setw(18) << corners[i](0)
                    << std::setw(18) << -corners[i](1) << "\n";
            }
            gpl << std::setw(18) << corners.front()(0)
                << std::setw(18) << -corners.front()(1) << "\n\n";
        }
    }

    elements_m.sort(ClassicField::SortAsc);
}

void OpalBeamline::save3DLattice() {
    if (Ippl::myNode() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    elements_m.sort([](const ClassicField& a, const ClassicField& b) {
            return a.order_m < b.order_m;
        });

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::ofstream pos;
    std::string fileName = "data/" + OpalData::getInstance()->getInputBasename() + "_ElementPositions.txt";
    if (Options::openMode == Options::APPEND && boost::filesystem::exists(fileName)) {
        pos.open(fileName, std::ios_base::app);
    } else {
        pos.open(fileName);
    }

    MeshGenerator mesh;
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        CoordinateSystemTrafo toBegin = element->getEdgeToBegin() * (*it).getCoordTransformationTo();
        CoordinateSystemTrafo toEnd = element->getEdgeToEnd() * (*it).getCoordTransformationTo();
        Vector_t entry3D = toBegin.getOrigin();
        Vector_t exit3D = toEnd.getOrigin();

        mesh.add(*(element.get()));

        if (element->getType() == ElementBase::SBEND ||
            element->getType() == ElementBase::RBEND) {

            Bend * bendElement = static_cast<Bend*>(element.get());
            std::vector<Vector_t> designPath = bendElement->getDesignPath();
            CoordinateSystemTrafo toEnd = bendElement->getBeginToEnd_local() * (*it).getCoordTransformationTo();
            Vector_t exit3D = toEnd.getOrigin();

            unsigned int size = designPath.size();
            unsigned int minNumSteps = std::max(20.0,
                                                std::abs(bendElement->getBendAngle() / Physics::pi * 180));
            unsigned int frequency = std::floor((double)size / minNumSteps);

            pos << std::setw(30) << std::left << std::string("\"ENTRY EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2)
                << std::setw(18) << std::setprecision(10) << entry3D(0)
                << std::setw(18) << std::setprecision(10) << entry3D(1)
                << "\n";

            Vector_t position = (*it).getCoordTransformationTo().transformFrom(designPath.front());
            pos << std::setw(30) << std::left << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2)
                << std::setw(18) << std::setprecision(10) << position(0)
                << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            for (unsigned int i = frequency; i + 1 < size; i += frequency) {

                Vector_t position = (*it).getCoordTransformationTo().transformFrom(designPath[i]);
                pos << std::setw(30) << std::left << std::string("\"MID: ") + element->getName() + std::string("\"")
                    << std::setw(18) << std::setprecision(10) << position(2)
                    << std::setw(18) << std::setprecision(10) << position(0)
                    << std::setw(18) << std::setprecision(10) << position(1)
                    << endl;
            }

            position = (*it).getCoordTransformationTo().transformFrom(designPath.back());
            pos << std::setw(30) << std::left << std::string("\"END: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2)
                << std::setw(18) << std::setprecision(10) << position(0)
                << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            pos << std::setw(30) << std::left << std::string("\"EXIT EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2)
                << std::setw(18) << std::setprecision(10) << exit3D(0)
                << std::setw(18) << std::setprecision(10) << exit3D(1)
                << std::endl;
        } else {
            pos << std::setw(30) << std::left << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2)
                << std::setw(18) << std::setprecision(10) << entry3D(0)
                << std::setw(18) << std::setprecision(10) << entry3D(1)
                << "\n";

            pos << std::setw(30) << std::left << std::string("\"END: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2)
                << std::setw(18) << std::setprecision(10) << exit3D(0)
                << std::setw(18) << std::setprecision(10) << exit3D(1)
                << std::endl;
        }
    }
    elements_m.sort(ClassicField::SortAsc);
    mesh.write(OpalData::getInstance()->getInputBasename());
}

namespace {
    std::string parseInput() {

        std::ifstream in(OpalData::getInstance()->getInputFn());
        std::string source("");
        std::string str;
        char testBit;
        const std::string commentFormat("");
        const boost::regex empty("^[ \t]*$");
        const boost::regex lineEnd(";");
        const std::string lineEndFormat(";\n");
        const boost::regex cppCommentExpr("//.*");
        const boost::regex cCommentExpr("/\\*.*?\\*/"); // "/\\*(?>[^*/]+|\\*[^/]|/[^*])*(?>(?R)(?>[^*/]+|\\*[^/]|/[^*])*)*\\*/"
        bool priorEmpty = true;

        in.get(testBit);
        while (!in.eof()) {
            in.putback(testBit);

            std::getline(in, str);
            str = boost::regex_replace(str, cppCommentExpr, commentFormat);
            str = boost::regex_replace(str, empty, commentFormat);
            if (str.size() > 0) {
                source += str;// + '\n';
                priorEmpty = false;
            } else if (!priorEmpty) {
                source += "##EMPTY_LINE##";
                priorEmpty = true;
            }

            in.get(testBit);
        }

        source = boost::regex_replace(source, cCommentExpr, commentFormat);
        source = boost::regex_replace(source, lineEnd, lineEndFormat, boost::match_default | boost::format_all);

        return source;
    }

    unsigned int getMinimalSignificantDigits(double num, const unsigned int maxDigits) {
        char buf[32];
        snprintf(buf, 32, "%.*f", maxDigits + 1, num);
        string numStr(buf);
        unsigned int length = numStr.length();

        unsigned int numDigits = maxDigits;
        unsigned int i = 2;
        while (i < maxDigits + 1 && numStr[length - i] == '0') {
            --numDigits;
            ++ i;
        }

        return numDigits;
    }

    std::string round2string(double num, const unsigned int maxDigits) {
        char buf[64];

        snprintf(buf, 64, "%.*f", getMinimalSignificantDigits(num, maxDigits), num);

        return std::string(buf);
    }
}

void OpalBeamline::save3DInput() {
    if (Ippl::myNode() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::string input = parseInput();
    std::ofstream pos("data/" + OpalData::getInstance()->getInputBasename() + "_3D.opal");

    for (; it != end; ++ it) {
        std::string element = (*it).getElement()->getName();
        const boost::regex replacePSI("(" + element + "\\s*:[^\\n]*)PSI\\s*=[^,;]*,?", boost::regex::icase);
        input = boost::regex_replace(input, replacePSI, "\\1\\2");

        const boost::regex replaceELEMEDGE("(" + element + "\\s*:[^\\n]*)ELEMEDGE\\s*=[^,;]*(.)", boost::regex::icase);

        CoordinateSystemTrafo cst = (*it).getCoordTransformationTo();
        Vector_t origin = cst.getOrigin();
        Vector_t orient = Util::getTaitBryantAngles(cst.getRotation().conjugate(), element);
        for (unsigned int d = 0; d < 3; ++ d)
            orient(d) *= Physics::rad2deg;

        std::string x = (std::abs(origin(0)) > 1e-10? "X = " + round2string(origin(0), 10) + ", ": "");
        std::string y = (std::abs(origin(1)) > 1e-10? "Y = " + round2string(origin(1), 10) + ", ": "");
        std::string z = (std::abs(origin(2)) > 1e-10? "Z = " + round2string(origin(2), 10) + ", ": "");

        std::string theta = (orient(0) > 1e-10? "THETA = " + round2string(orient(0), 6) + " * PI / 180, ": "");
        std::string phi = (orient(1) > 1e-10? "PHI = " + round2string(orient(1), 6) + " * PI / 180, ": "");
        std::string psi = (orient(2) > 1e-10? "PSI = " + round2string(orient(2), 6) + " * PI / 180, ": "");
        std::string coordTrafo = x + y + z + theta + phi + psi;
        if (coordTrafo.length() > 2) {
            coordTrafo = coordTrafo.substr(0, coordTrafo.length() - 2); // remove last ', '
        }

        std::string position = ("\\1" + coordTrafo + "\\2");

        input = boost::regex_replace(input, replaceELEMEDGE, position);

        if ((*it).getElement()->getType() == ElementBase::RBEND ||
            (*it).getElement()->getType() == ElementBase::SBEND) {
            const Bend* dipole = static_cast<const Bend*>((*it).getElement().get());
            double angle = dipole->getBendAngle();
            double E1 = dipole->getEntranceAngle();
            double E2 = dipole->getExitAngle();

            const boost::regex angleR("(" + element + "\\s*:[^\\n]*ANGLE\\s*=)[^,;]*(.)");
            const std::string angleF("\\1 " + round2string(angle * 180 / Physics::pi, 6) + " / 180 * PI\\2");
            const boost::regex E1R("(" + element + "\\s*:[^\\n]*E1\\s*=)[^,;]*(.)");
            const std::string E1F("\\1 " + round2string(E1 * 180 / Physics::pi, 6) + " / 180 * PI\\2");
            const boost::regex E2R("(" + element + "\\s*:[^\\n]*E2\\s*=)[^,;]*(.)");
            const std::string E2F("\\1 " + round2string(E2 * 180 / Physics::pi, 6) + " / 180 * PI\\2");
            const boost::regex noRotation("(" + element + "\\s*:[^\\n]*),\\s*ROTATION\\s*=[^,;]*(.)");
            const std::string noRotationFormat("\\1\\2  ");

            input = boost::regex_replace(input, angleR, angleF);
            input = boost::regex_replace(input, E1R, E1F);
            input = boost::regex_replace(input, E2R, E2F);
            input = boost::regex_replace(input, noRotation, noRotationFormat);
        }
    }

    const boost::regex empty("##EMPTY_LINE##");
    const std::string emptyFormat("\n");
    input = boost::regex_replace(input, empty, emptyFormat);

    pos << input << std::endl;
}

void OpalBeamline::activateElements() {
    auto it = elements_m.begin();
    const auto end = elements_m.end();

    double designEnergy = 0.0;
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        if (element->getType() == ElementBase::SBEND ||
            element->getType() == ElementBase::RBEND) {
            Bend * bendElement = static_cast<Bend*>(element.get());
            designEnergy = bendElement->getDesignEnergy() * 1e-6;
        }
        (*it).setOn(designEnergy);
        // element->goOnline(designEnergy);
    }
}