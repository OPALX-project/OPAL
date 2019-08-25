#include "Structure/MeshGenerator.h"
#include "Physics/Physics.h"
#include "Utilities/Util.h"
#include "AbsBeamline/Bend2D.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/Multipole.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

MeshGenerator::MeshGenerator():
    elements_m()
{ }

void MeshGenerator::add(const ElementBase &element) {
    double start = 0.0;

    MeshData mesh;
    if (element.getType() == ElementBase::SBEND ||
        element.getType() == ElementBase::RBEND) {

        const Bend2D* dipole = static_cast<const Bend2D*>(&element);
        mesh = dipole->getSurfaceMesh();
        mesh.type_m = DIPOLE;
    } else if (element.getType() == ElementBase::RBEND3D) {
        const RBend3D* dipole = static_cast<const RBend3D*>(&element);
        mesh = dipole->getSurfaceMesh();
        mesh.type_m = DIPOLE;
    } else if (element.getType() == ElementBase::DRIFT) {
        return;
    } else {
        double end, length;
        element.getElementDimensions(start, end);
        length = end - start;
        auto apert = element.getAperture();

        switch (apert.first) {
        case ElementBase::RECTANGULAR:
        case ElementBase::CONIC_RECTANGULAR:
            mesh = getBox(length, apert.second[0], apert.second[1], apert.second[2]);
            break;
        case ElementBase::ELLIPTICAL:
        case ElementBase::CONIC_ELLIPTICAL:
            mesh = getCylinder(length, apert.second[0], apert.second[1], apert.second[2]);
            break;
        default:
            return;
        }

        switch(element.getType()) {
        case ElementBase::MULTIPOLE:
            switch(static_cast<const Multipole*>(&element)->getMaxNormalComponentIndex()) {
            case 1:
                mesh.type_m = DIPOLE;
                break;
            case 2:
                mesh.type_m = QUADRUPOLE;
                break;
            case 3:
                mesh.type_m = SEXTUPOLE;
                break;
            case 4:
                mesh.type_m = OCTUPOLE;
                break;
            default:
                break;
            }
            break;
        case ElementBase::RFCAVITY:
        case ElementBase::VARIABLERFCAVITY:
        case ElementBase::TRAVELINGWAVE:
            mesh.type_m = RFCAVITY;
            break;
        case ElementBase::SOLENOID:
            mesh.type_m = SOLENOID;
            break;
        default:
            mesh.type_m = OTHER;
        }
    }

    CoordinateSystemTrafo trafo = element.getCSTrafoGlobal2Local().inverted();
    Vector_t z = trafo.rotateTo(Vector_t(0, 0, 1));
    for (unsigned int i = 0; i < mesh.vertices_m.size(); ++ i) {
        mesh.vertices_m[i] = trafo.transformTo(mesh.vertices_m[i]) + start * z;
    }

    for (unsigned int i = 0; i < mesh.decorations_m.size(); ++ i) {
        mesh.decorations_m[i].first = trafo.transformTo(mesh.decorations_m[i].first) + start * z;
        mesh.decorations_m[i].second = trafo.transformTo(mesh.decorations_m[i].second) + start * z;
    }

    elements_m.push_back(mesh);
}

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>

void MeshGenerator::write(const std::string &fname) {
    std::ofstream out("data/" + fname + "_ElementPositions.py");
    const char *buffer;
    const std::string indent("    ");

    out << std::fixed << std::setprecision(6);

    out << "import os, sys, argparse, math, base64, zlib, struct\n\n";
    out << "if sys.version_info < (3,0):\n";
    out << "    range = xrange\n\n";

    std::stringstream vertices_ascii;
    std::ostringstream vertices_compressed;
    out << "vertices = []\n";
    out << "numVertices = [";
    for (auto &element: elements_m) {
        auto &vertices = element.vertices_m;
        out << vertices.size() << ", ";
        for (unsigned int i = 0; i < vertices.size(); ++ i) {
            buffer = reinterpret_cast<const char*>(&(vertices[i](0)));
            vertices_ascii.write(buffer, 3 * sizeof(double));
        }
    }
    out.seekp(-2, std::ios_base::end);
    out << "]\n";

    {
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        in.push(boost::iostreams::zlib_compressor());
        in.push(vertices_ascii);
        boost::iostreams::copy(in, vertices_compressed);
    }

    std::string vertices_base64 = Util::base64_encode(vertices_compressed.str());
    out << "vertices_base64 = \"" << vertices_base64 << "\"\n";
    out << "triangles = [";
    for (auto &element: elements_m) {
        out << "[ ";
        auto &triangles = element.triangles_m;
        for (unsigned int i = 0; i < triangles.size(); ++ i) {
            out << triangles[i](0) << ", " << triangles[i](1) << ", " << triangles[i](2) << ", ";
        }
        out.seekp(-2, std::ios_base::end);
        out << "], ";
    }
    out.seekp(-2, std::ios_base::end);
    out << "]\n";


    out << "decoration = [";
    for (auto &element: elements_m) {
        out << "[ ";
        for (auto &decoration: element.decorations_m) {
            out << (decoration.first)(0) << ", " << (decoration.first)(1) << ", " << (decoration.first)(2) << ", "
                << (decoration.second)(0) << ", " << (decoration.second)(1) << ", " << (decoration.second)(2) << ", ";
        }
        if (element.decorations_m.size() > 0)
            out.seekp(-2, std::ios_base::end);
        out << "], ";
    }
    if (elements_m.size() > 0)
        out.seekp(-2, std::ios_base::end);
    out << "]\n\n";

    out << "color = [";
    for (auto &element: elements_m) {
        out << element.type_m << ", ";
    }
    out.seekp(-2, std::ios_base::end);
    out << "]\n\n";

    std::stringstream index_ascii;
    std::ostringstream index_compressed;
    index_ascii << "<!DOCTYPE html>\n"
                << "<html>\n"
                << "    <head>\n"
                << "        <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n"
                << "\n"
                << "        <title>Babylon.js sample code</title>\n"
                << "        <!-- Babylon.js -->\n"
                << "        <script src=\"http://cdn.babylonjs.com/hand.minified-1.2.js\"></script>\n"
                << "        <script src=\"http://cdn.babylonjs.com/cannon.js\"></script>\n"
                << "        <script src=\"http://cdn.babylonjs.com/oimo.js\"></script>\n"
                << "        <script src=\"http://cdn.babylonjs.com/babylon.js\"></script>\n"
                << "        <style>\n"
                << "            html, body {\n"
                << "                overflow: hidden;\n"
                << "                width: 100%;\n"
                << "                height: 100%;\n"
                << "                margin: 0;\n"
                << "                padding: 0;\n"
                << "            }\n"
                << "\n"
                << "            #renderCanvas {\n"
                << "                width: 100%;\n"
                << "                height: 100%;\n"
                << "                touch-action: none;\n"
                << "            }\n"
                << "        </style>\n"
                << "    </head>\n"
                << "<body>\n"
                << "    <canvas id=\"renderCanvas\"></canvas>\n"
                << "    <script>\n"
                << "        var canvas = document.getElementById(\"renderCanvas\");\n"
                << "        var engine = new BABYLON.Engine(canvas, true);\n"
                << "\n"
                << "        var createScene = function () {\n"
                << "            var scene = new BABYLON.Scene(engine);\n"
                << "\n"
                << "            //Adding a light\n"
                << "                var light = new BABYLON.PointLight(\"Omni\", new BABYLON.Vector3(0, -100, 0), scene);\n"
                << "                //light.diffuse = new BABYLON.Color3(0.8, 0.8, 0.8);\n"
                << "                light.specular = new BABYLON.Color3(0.1, 0.1, 0.1);\n"
                << "                //light.groundColor = new BABYLON.Color3(0, 0, 0);\n"
                << "\n"
                << "            //Adding an Arc Rotate Camera\n"
                << "                var camera = new BABYLON.ArcRotateCamera(\"Camera\", 0.0, Math.PI, 50, BABYLON.Vector3.Zero(), scene);\n"
                << "            camera.attachControl(canvas, false);\n"
                << "                camera.wheelPrecision = 20;\n"
                << "\n"
                << "            var mymesh = ##DATA##;\n"
                << "            // The first parameter can be used to specify which mesh to import. Here we import all meshes\n"
                << "            BABYLON.SceneLoader.ImportMesh(\"\", \"\", mymesh, scene, function (newMeshes) {\n"
                << "                // Set the target of the camera to the first imported mesh\n"
                << "                camera.target = newMeshes[0];\n"
                << "\n"
                << "                        var material1 = new BABYLON.StandardMaterial(\"texture1\", scene);\n"
                << "                        //material1.wireframe = true;\n"
                << "                newMeshes[0].material = material1;\n"
                << "                newMeshes[0].material.emissiveColor = new BABYLON.Color3(0.2, 0.2, 0.2);\n"
                << "            });\n"
                << "\n"
                << "            return scene;\n"
                << "        }\n"
                << "\n"
                << "\n"
                << "        var scene = createScene();\n"
                << "\n"
                << "        engine.runRenderLoop(function () {\n"
                << "            scene.render();\n"
                << "        });\n"
                << "\n"
                << "        // Resize\n"
                << "        window.addEventListener(\"resize\", function () {\n"
                << "            engine.resize();\n"
                << "        });\n"
                << "    </script>\n"
                << "</body>\n"
                << "</html>";
    {
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        in.push(boost::iostreams::zlib_compressor());
        in.push(index_ascii);
        boost::iostreams::copy(in, index_compressed);
    }

    out << "index_base64 = '" << Util::base64_encode(index_compressed.str()) << "'\n\n";

    out << "def decodeVertices():\n";
    out << indent << "vertices_binary = zlib.decompress(base64.b64decode(vertices_base64))\n";
    out << indent << "k = 0\n";
    out << indent << "for i in range(len(numVertices)):\n";
    out << indent << indent << "current = []\n";
    out << indent << indent << "for j in range(numVertices[i] * 3):\n";
    out << indent << indent << indent << "current.append(float(struct.unpack('=d', vertices_binary[k:k+8])[0]))\n";
    out << indent << indent << indent << "k += 8\n";
    out << indent << indent << "vertices.append(current)\n\n";

    out << "def dot(a, b):\n";
    out << indent << "return sum(a[i]*b[i] for i in range(len(a)))\n\n";

    out << "def normalize(a):\n";
    out << indent << "length = math.sqrt(dot(a, a))\n";
    out << indent << "for i in range(3):\n";
    out << indent << indent << "a[i]/=length\n";
    out << indent << "return a\n\n";

    out << "def distance(a, b):\n";
    out << indent << "c = [0] * 2\n";
    out << indent << "for i in range(len(a)):\n";
    out << indent << indent << "c[i] = a[i] - b[i]\n";
    out << indent << "return math.sqrt(dot(c, c))\n\n";

    out << "def cross(a, b):\n";
    out << indent << "c = [0, 0, 0]\n";
    out << indent << "c[0] = a[1]*b[2] - a[2]*b[1]\n";
    out << indent << "c[1] = a[2]*b[0] - a[0]*b[2]\n";
    out << indent << "c[2] = a[0]*b[1] - a[1]*b[0]\n";
    out << indent << "return c\n\n";

    out << "class Quaternion:\n";
    out << indent << "def __init__(self, values):\n";
    out << indent << indent << "self.R = values\n\n";

    out << indent << "def __str__(self):\n";
    out << indent << indent << "return str(self.R)\n\n";

    out << indent << "def __mul__(self, other):\n";
    out << indent << indent << "imagThis = self.imag()\n";
    out << indent << indent << "imagOther = other.imag()\n";
    out << indent << indent << "cro = cross(imagThis, imagOther)\n\n";

    out << indent << indent << "ret = ([self.R[0] * other.R[0] - dot(imagThis, imagOther)] +\n";
    out << indent << indent << indent << "   [self.R[0] * imagOther[0] + other.R[0] * imagThis[0] + cro[0],\n";
    out << indent << indent << indent << indent << "self.R[0] * imagOther[1] + other.R[0] * imagThis[1] + cro[1],\n";
    out << indent << indent << indent << indent << "self.R[0] * imagOther[2] + other.R[0] * imagThis[2] + cro[2]])\n\n";

    out << indent << indent << "return Quaternion(ret)\n\n";

    out << indent << "def imag(self):\n";
    out << indent << indent << "return self.R[1:]\n\n";

    out << indent << "def conjugate(self):\n";
    out << indent << indent << "ret = [0] * 4\n";
    out << indent << indent << "ret[0] = self.R[0]\n";
    out << indent << indent << "ret[1] = -self.R[1]\n";
    out << indent << indent << "ret[2] = -self.R[2]\n";
    out << indent << indent << "ret[3] = -self.R[3]\n\n";

    out << indent << indent << "return Quaternion(ret)\n\n";

    out << indent << "def inverse(self):\n";
    out << indent << indent << "return self.conjugate()\n\n";

    out << indent << "def rotate(self, vec3D):\n";
    out << indent << indent << "quat = Quaternion([0.0] + vec3D)\n";
    out << indent << indent << "conj = self.conjugate()\n\n";

    out << indent << indent << "ret = self * (quat * conj)\n";
    out << indent << indent << "return ret.R[1:]\n\n";

    out << "def getQuaternion(u, ref):\n";
    out << indent << "u = normalize(u)\n";
    out << indent << "ref = normalize(ref)\n";
    out << indent << "axis = cross(u, ref)\n";
    out << indent << "normAxis = math.sqrt(dot(axis, axis))\n\n";

    out << indent << "if normAxis < 1e-12:\n";
    out << indent << indent << "if math.fabs(dot(u, ref) - 1.0) < 1e-12:\n";
    out << indent << indent << indent << "return Quaternion([1.0, 0.0, 0.0, 0.0])\n\n";

    out << indent << indent << "return Quaternion([0.0, 1.0, 0.0, 0.0])\n\n";

    out << indent << "axis = normalize(axis)\n";
    out << indent << "cosAngle = math.sqrt(0.5 * (1 + dot(u, ref)))\n";
    out << indent << "sinAngle = math.sqrt(1 - cosAngle * cosAngle)\n\n";

    out << indent << "return Quaternion([cosAngle, sinAngle * axis[0], sinAngle * axis[1], sinAngle * axis[2]])\n\n";

    out << "def exportVTK():\n";
    out << indent << "vertices_str = \"\"\n";
    out << indent << "triangles_str = \"\"\n";
    out << indent << "color_str = \"\"\n";
    out << indent << "cellTypes_str = \"\"\n";
    out << indent << "startIdx = 0\n";
    out << indent << "vertexCounter = 0\n";
    out << indent << "cellCounter = 0\n";
    out << indent << "lookup_table = []\n";
    out << indent << "lookup_table.append([0.5, 0.5, 0.5, 0.5])\n";
    out << indent << "lookup_table.append([1.0, 0.847, 0.0, 1.0])\n";
    out << indent << "lookup_table.append([1.0, 0.0, 0.0, 1.0])\n";
    out << indent << "lookup_table.append([0.537, 0.745, 0.525, 1.0])\n";
    out << indent << "lookup_table.append([0.5, 0.5, 0.0, 1.0])\n";
    out << indent << "lookup_table.append([1.0, 0.541, 0.0, 1.0])\n";
    out << indent << "lookup_table.append([0.0, 0.0, 1.0, 1.0])\n\n";

    out << indent << "decodeVertices()\n\n";

    out << indent << "for i in range(len(vertices)):\n";
    out << indent << indent << "for j in range(0, len(vertices[i]), 3):\n";
    out << indent << indent << indent << "vertices_str += (\"%f %f %f\\n\" %(vertices[i][j], vertices[i][j+1], vertices[i][j+2]))\n";
    out << indent << indent << indent << "vertexCounter += 1\n\n";

    out << indent << indent << "for j in range(0, len(triangles[i]), 3):\n";
    out << indent << indent << indent << "triangles_str += (\"3 %d %d %d\\n\" % (triangles[i][j] + startIdx, triangles[i][j+1] + startIdx, triangles[i][j+2] + startIdx))\n";
    out << indent << indent << indent << "cellTypes_str += \"5\\n\"\n";
    out << indent << indent << indent << "tmp_color = lookup_table[color[i]]\n";
    out << indent << indent << indent << "color_str += (\"%f %f %f %f\\n\" % (tmp_color[0], tmp_color[1], tmp_color[2], tmp_color[3]))\n";
    out << indent << indent << indent << "cellCounter += 1\n";

    out << indent << indent << "startIdx = vertexCounter\n\n";

    out << indent << "fh = open('" << fname << "_ElementPositions.vtk','w')\n";
    out << indent << "fh.write(\"# vtk DataFile Version 2.0\\n\")\n";
    out << indent << "fh.write(\"test\\nASCII\\n\\n\")\n";
    out << indent << "fh.write(\"DATASET UNSTRUCTURED_GRID\\n\")\n";
    out << indent << "fh.write(\"POINTS \" + str(vertexCounter) + \" float\\n\")\n";
    out << indent << "fh.write(vertices_str)\n";
    out << indent << "fh.write(\"CELLS \" + str(cellCounter) + \" \" + str(cellCounter * 4) + \"\\n\")\n";
    out << indent << "fh.write(triangles_str)\n";
    out << indent << "fh.write(\"CELL_TYPES \" + str(cellCounter) + \"\\n\")\n";
    out << indent << "fh.write(cellTypes_str)\n";
    out << indent << "fh.write(\"CELL_DATA \" + str(cellCounter) + \"\\n\")\n";
    out << indent << "fh.write(\"COLOR_SCALARS type 4\\n\")\n";
    out << indent << "fh.write(color_str + \"\\n\")\n";
    out << indent << "fh.close()\n\n";

    out << "def getNormal(tri_vertices):\n";
    out << indent << "vec1 = [tri_vertices[1][0] - tri_vertices[0][0],\n";
    out << indent << "        tri_vertices[1][1] - tri_vertices[0][1],\n";
    out << indent << "        tri_vertices[1][2] - tri_vertices[0][2]]\n";
    out << indent << "vec2 = [tri_vertices[2][0] - tri_vertices[0][0],\n";
    out << indent << "        tri_vertices[2][1] - tri_vertices[0][1],\n";
    out << indent << "        tri_vertices[2][2] - tri_vertices[0][2]]\n";
    out << indent << "return normalize(cross(vec1,vec2))\n\n";

    out << "def exportWeb(bgcolor):\n";
    // out << indent << "if not os.path.exists('scenes'):\n";
    // out << indent << indent << "os.makedirs('scenes')\n";
    // out << indent << "fh = open('scenes/" << fname << "_ElementPositions.babylon','w')\n";
    // out << indent << "indent = \"\";\n\n";

    out << indent << "lookup_table = []\n";
    out << indent << "lookup_table.append([0.5, 0.5, 0.5])\n";
    out << indent << "lookup_table.append([1.0, 0.847, 0.0])\n";
    out << indent << "lookup_table.append([1.0, 0.0, 0.0])\n";
    out << indent << "lookup_table.append([0.537, 0.745, 0.525])\n";
    out << indent << "lookup_table.append([0.5, 0.5, 0.0])\n";
    out << indent << "lookup_table.append([1.0, 0.541, 0.0])\n";
    out << indent << "lookup_table.append([0.0, 0.0, 1.0])\n\n";

    out << indent << "decodeVertices()\n\n";

    out << indent << "mesh = \"'data:\"\n";
    out << indent << "mesh += \"{\"\n";
    out << indent << "mesh += \"\\\"autoClear\\\":true,\"\n";
    out << indent << "mesh += \"\\\"clearColor\\\":[0.0000,0.0000,0.0000],\"\n";
    out << indent << "mesh += \"\\\"ambientColor\\\":[0.0000,0.0000,0.0000],\"\n";
    out << indent << "mesh += \"\\\"gravity\\\":[0.0000,-9.8100,0.0000],\"\n";
    out << indent << "mesh += \"\\\"cameras\\\":[\"\n";
    out << indent << "mesh += \"{\"\n";
    out << indent << "mesh += \"\\\"name\\\":\\\"Camera\\\",\"\n";
    out << indent << "mesh += \"\\\"id\\\":\\\"Camera\\\",\"\n";
    out << indent << "mesh += \"\\\"position\\\":[21.7936,2.2312,-85.7292],\"\n";
    out << indent << "mesh += \"\\\"rotation\\\":[0.0432,-0.1766,-0.0668],\"\n";
    out << indent << "mesh += \"\\\"fov\\\":0.8578,\"\n";
    out << indent << "mesh += \"\\\"minZ\\\":10.0000,\"\n";
    out << indent << "mesh += \"\\\"maxZ\\\":10000.0000,\"\n";
    out << indent << "mesh += \"\\\"speed\\\":1.0000,\"\n";
    out << indent << "mesh += \"\\\"inertia\\\":0.9000,\"\n";
    out << indent << "mesh += \"\\\"checkCollisions\\\":false,\"\n";
    out << indent << "mesh += \"\\\"applyGravity\\\":false,\"\n";
    out << indent << "mesh += \"\\\"ellipsoid\\\":[0.2000,0.9000,0.2000]\"\n";
    out << indent << "mesh += \"}],\"\n";
    out << indent << "mesh += \"\\\"activeCamera\\\":\\\"Camera\\\",\"\n";
    out << indent << "mesh += \"\\\"lights\\\":[\"\n";
    out << indent << "mesh += \"{\"\n";
    out << indent << "mesh += \"\\\"name\\\":\\\"Lamp\\\",\"\n";
    out << indent << "mesh += \"\\\"id\\\":\\\"Lamp\\\",\"\n";
    out << indent << "mesh += \"\\\"type\\\":0.0000,\"\n";
    out << indent << "mesh += \"\\\"position\\\":[4.0762,34.9321,-63.5788],\"\n";
    out << indent << "mesh += \"\\\"intensity\\\":1.0000,\"\n";
    out << indent << "mesh += \"\\\"diffuse\\\":[1.0000,1.0000,1.0000],\"\n";
    out << indent << "mesh += \"\\\"specular\\\":[1.0000,1.0000,1.0000]\"\n";
    out << indent << "mesh += \"}],\"\n";
    out << indent << "mesh += \"\\\"materials\\\":[],\"\n";
    out << indent << "mesh += \"\\\"meshes\\\":[\"\n\n";

    out << indent << "for i in range(len(triangles)):\n";
    out << indent << indent << "vertex_list = []\n";
    out << indent << indent << "indices_list = []\n";
    out << indent << indent << "normals_list = []\n";
    out << indent << indent << "color_list = []\n";
    out << indent << indent << "for j in range(0, len(triangles[i]), 3):\n";
    out << indent << indent << indent << "tri_vertices = []\n";
    out << indent << indent << indent << "idcs = triangles[i][j:j + 3]\n";
    out << indent << indent << indent << "tri_vertices.append(vertices[i][3 * idcs[0]:3 * (idcs[0] + 1)])\n";
    out << indent << indent << indent << "tri_vertices.append(vertices[i][3 * idcs[1]:3 * (idcs[1] + 1)])\n";
    out << indent << indent << indent << "tri_vertices.append(vertices[i][3 * idcs[2]:3 * (idcs[2] + 1)])\n";
    out << indent << indent << indent << "indices_list.append(','.join(str(n) for n in range(len(vertex_list),len(vertex_list) + 3)))\n";
    out << indent << indent << indent << "# left hand order!\n";
    out << indent << indent << indent << "vertex_list.append(','.join(\"%.5f\" % (round(n,5) + 0.0) for n in tri_vertices[0]))\n";
    out << indent << indent << indent << "vertex_list.append(','.join(\"%.5f\" % (round(n,5) + 0.0) for n in tri_vertices[2]))\n";
    out << indent << indent << indent << "vertex_list.append(','.join(\"%.5f\" % (round(n,5) + 0.0) for n in tri_vertices[1]))\n";
    out << indent << indent << indent << "normal = getNormal(tri_vertices)\n";
    out << indent << indent << indent << "normals_list.append(','.join(\"%.5f\" % (round(n,5) + 0.0) for n in normal * 3))\n";
    out << indent << indent << indent << "color_list.append(','.join([str(n) for n in lookup_table[color[i]]] * 3))\n";
    out << indent << indent << "mesh += \"{\"\n";
    out << indent << indent << "mesh += \"\\\"name\\\":\\\"element_\" + str(i) + \"\\\",\"\n";
    out << indent << indent << "mesh += \"\\\"id\\\":\\\"element_\" + str(i) + \"\\\",\"\n";
    out << indent << indent << "mesh += \"\\\"position\\\":[0.0,0.0,0.0],\"\n";
    out << indent << indent << "mesh += \"\\\"rotation\\\":[0.0,0.0,0.0],\"\n";
    out << indent << indent << "mesh += \"\\\"scaling\\\":[1.0,1.0,1.0],\"\n";
    out << indent << indent << "mesh += \"\\\"isVisible\\\":true,\"\n";
    out << indent << indent << "mesh += \"\\\"isEnabled\\\":true,\"\n";
    out << indent << indent << "mesh += \"\\\"useFlatShading\\\":false,\"\n";
    out << indent << indent << "mesh += \"\\\"checkCollisions\\\":false,\"\n";
    out << indent << indent << "mesh += \"\\\"billboardMode\\\":0,\"\n";
    out << indent << indent << "mesh += \"\\\"receiveShadows\\\":false,\"\n";
    out << indent << indent << "mesh += \"\\\"positions\\\":[\" + ','.join(vertex_list) + \"],\"\n";
    out << indent << indent << "mesh += \"\\\"normals\\\":[\" + ','.join(normals_list) + \"],\"\n";
    out << indent << indent << "mesh += \"\\\"indices\\\":[\" + ','.join(indices_list) + \"],\"\n";
    out << indent << indent << "mesh += \"\\\"colors\\\":[\" + ','.join(color_list) + \"],\"\n";
    out << indent << indent << "mesh += \"\\\"subMeshes\\\":[\"\n";
    out << indent << indent << "mesh += \"{\"\n";
    out << indent << indent << "mesh += \"\\\"materialIndex\\\":0,\"\n";
    out << indent << indent << "mesh += \" \\\"verticesStart\\\":0,\"\n";
    out << indent << indent << "mesh += \" \\\"verticesCount\\\":\" + str(len(triangles[i])) + \",\"\n";
    out << indent << indent << "mesh += \" \\\"indexStart\\\":0,\"\n";
    out << indent << indent << "mesh += \" \\\"indexCount\\\":\" + str(len(triangles[i])) + \"\"\n";
    out << indent << indent << "mesh += \"}]\"\n";
    out << indent << indent << "mesh += \"}\"\n";
    out << indent << indent << "mesh += \",\"\n\n";

    out << indent << indent << "del normals_list[:]\n";
    out << indent << indent << "del vertex_list[:]\n";
    out << indent << indent << "del color_list[:]\n";
    out << indent << indent << "del indices_list[:]\n\n";

    out << indent << "mesh = mesh[:-1] + \"]\"\n";
    out << indent << "mesh += \"}'\"\n";

    out << indent << "index_compressed = base64.b64decode(index_base64)\n";
    out << indent << "index = str(zlib.decompress(index_compressed))\n";
    out << indent << "if (len(bgcolor) == 3):\n";
    out << indent << indent << "mesh += \";\\n            \"\n";
    out << indent << indent << "mesh += \"scene.clearColor = new BABYLON.Color3(%f, %f, %f)\" % (bgcolor[0], bgcolor[1], bgcolor[2])\n\n";
    out << indent << "index = index.replace('##DATA##', mesh)\n";
    out << indent << "fh = open('" << fname << "_ElementPositions.html','w')\n";
    out << indent << "fh.write(index)\n";
    out << indent << "fh.close()\n\n";

    out << "def computeMinAngle(idx, curAngle, positions, connections, check):\n";
    out << indent << "matrix = [[-math.cos(curAngle), -math.sin(curAngle)], [math.sin(curAngle), -math.cos(curAngle)]]\n\n";

    out << indent << "minAngle = 2 * math.pi\n";
    out << indent << "nextIdx = -1\n\n";

    out << indent << "for j in connections[idx]:\n";
    out << indent << indent << "direction = [positions[j][0] - positions[idx][0],\n";
    out << indent << indent << indent << "         positions[j][1] - positions[idx][1]]\n";
    out << indent << indent << "direction = [dot(matrix[0],direction), dot(matrix[1],direction)]\n\n";

    out << indent << indent << "if math.fabs(dot([1.0, 0.0], direction) / distance(positions[j], positions[idx]) - 1.0) < 1e-6: continue\n\n";

    out << indent << indent << "angle = math.fmod(math.atan2(direction[1],direction[0]) + 2 * math.pi, 2 * math.pi)\n\n";

    out << indent << indent << "if angle < minAngle:\n";
    out << indent << indent << indent << "nextIdx = j\n";
    out << indent << indent << indent << "minAngle = angle\n";
    out << indent << indent << "if angle == minAngle and check:\n";
    out << indent << indent << indent << "dire =  [positions[j][0] - positions[idx][0],\n";
    out << indent << indent << indent << "         positions[j][1] - positions[idx][1]]\n";
    out << indent << indent << indent << "minA0 = math.atan2(dire[1], dire[0])\n";
    out << indent << indent << indent << "minA1 = computeMinAngle(nextIdx, minA0, positions, connections, False)\n";
    out << indent << indent << indent << "minA2 = computeMinAngle(j, minA0, positions, connections, False)\n";
    out << indent << indent << indent << "if minA2[1] < minA2[1]:\n";
    out << indent << indent << indent << indent << "nextIdx = j\n\n";

    out << indent << "if nextIdx == -1:\n";
    out << indent << indent << "nextIdx = connections[idx][0]\n\n";

    out << indent << "return (nextIdx, minAngle)\n\n";

    out << "def squashVertices(positionDict, connections):\n";
    out << indent << "removedItems = []\n";
    out << indent << "indices = [int(k) for k in positionDict.keys()]\n";
    out << indent << "idxChanges = indices\n";
    out << indent << "for i in indices:\n";
    out << indent << indent << "if i in removedItems:\n";
    out << indent << indent << indent << "continue\n";
    out << indent << indent << "for j in indices:\n";
    out << indent << indent << indent << "if j in removedItems or j == i:\n";
    out << indent << indent << indent << indent << "continue\n\n";

    out << indent << indent << indent << "if distance(positionDict[i], positionDict[j]) < 1e-6:\n";
    out << indent << indent << indent << indent << "connections[j] = list(set(connections[j]))\n";
    out << indent << indent << indent << indent << "if i in connections[j]:\n";
    out << indent << indent << indent << indent << indent << "connections[j].remove(i)\n";
    out << indent << indent << indent << indent << "if j in connections[i]:\n";
    out << indent << indent << indent << indent << indent << "connections[i].remove(j)\n\n";

    out << indent << indent << indent << indent << "connections[i].extend(connections[j])\n";
    out << indent << indent << indent << indent << "connections[i] = list(set(connections[i]))\n";
    out << indent << indent << indent << indent << "connections[i].sort()\n\n";

    out << indent << indent << indent << indent << "for k in connections.keys():\n";
    out << indent << indent << indent << indent << indent << "if k == i: continue\n";
    out << indent << indent << indent << indent << indent << "if j in connections[k]:\n";
    out << indent << indent << indent << indent << indent << indent << "idx = connections[k].index(j)\n";
    out << indent << indent << indent << indent << indent << indent << "connections[k][idx] = i\n\n";

    out << indent << indent << indent << indent << "idxChanges[j] = i\n";
    out << indent << indent << indent << indent << "removedItems.append(j)\n\n";

    out << indent << "for i in removedItems:\n";
    out << indent << indent << "del positionDict[i]\n";
    out << indent << indent << "del connections[i]\n\n";

    out << indent << "for i in connections.keys():\n";
    out << indent << indent << "connections[i] = list(set(connections[i]))\n";
    out << indent << indent << "connections[i].sort()\n\n";

    out << indent << "return idxChanges\n\n";

    out << "def computeLineEquations(positions, connections):\n";
    out << indent << "cons = set()\n";
    out << indent << "for i in connections.keys():\n";
    out << indent << indent << "for j in connections[i]:\n";
    out << indent << indent << indent << "cons.add((min(i, j), max(i, j)))\n\n";

    out << indent << "lineEquations = {}\n";
    out << indent << "for item in cons:\n";
    out << indent << indent << "a = (positions[item[1]][1] - positions[item[0]][1])\n";
    out << indent << indent << "b = -(positions[item[1]][0] - positions[item[0]][0])\n\n";

    out << indent << indent << "xm = 0.5 * (positions[item[0]][0] + positions[item[1]][0])\n";
    out << indent << indent << "ym = 0.5 * (positions[item[0]][1] + positions[item[1]][1])\n";
    out << indent << indent << "c = -(a * xm +  b * ym)\n\n";

    out << indent << indent << "lineEquations[item] = (a, b, c)\n\n";

    out << indent << "return lineEquations\n\n";

    out << "def checkPossibleSegmentIntersection(segment, positions, lineEquations, minAngle, lastIntersection):\n";
    out << indent << "item1 = (min(segment), max(segment))\n\n";

    out << indent << "(a1, b1, c1) = (0,0,0)\n";
    out << indent << "A = [0]*2\n";
    out << indent << "B = A\n\n";

    out << indent << "if segment[0] == None:\n";
    out << indent << indent << "A = lastIntersection\n";
    out << indent << indent << "B = positions[segment[1]]\n\n";

    out << indent << indent << "a1 = B[1] - A[1]\n";
    out << indent << indent << "b1 = -(B[0] - A[0])\n";
    out << indent << indent << "xm = 0.5 * (A[0] + B[0])\n";
    out << indent << indent << "ym = 0.5 * (A[1] + B[1])\n";
    out << indent << indent << "c1 = -(a1 * xm + b1 * ym)\n\n";

    out << indent << "else:\n";
    out << indent << indent << "A = positions[segment[0]]\n";
    out << indent << indent << "B = positions[segment[1]]\n\n";

    out << indent << indent << "(a1, b1, c1) = lineEquations[item1]\n\n";

    out << indent << indent << "if segment[1] < segment[0]:\n";
    out << indent << indent << indent << "(a1, b1, c1) = (-a1, -b1, -c1)\n\n";

    out << indent << "curAngle = math.atan2(a1, -b1)\n";
    out << indent << "matrix = [[-math.cos(curAngle), -math.sin(curAngle)], [math.sin(curAngle), -math.cos(curAngle)]]\n\n";

    out << indent << "origMinAngle = minAngle\n\n";

    out << indent << "segment1 = [B[0] - A[0], B[1] - A[1], 0.0]\n\n";

    out << indent << "intersectingSegments = []\n";
    out << indent << "distanceAB = distance(A, B)\n";
    out << indent << "for item2 in lineEquations.keys():\n";
    out << indent << indent << "item = item2\n";
    out << indent << indent << "C = positions[item[0]]\n";
    out << indent << indent << "D = positions[item[1]]\n\n";

    out << indent << indent << "if segment[0] == None:\n";
    out << indent << indent << indent << "if (segment[1] == item[0] or\n";
    out << indent << indent << indent << indent << "segment[1] == item[1]): continue\n";
    out << indent << indent << "else:\n";
    out << indent << indent << indent << "if (item1[0] == item[0] or\n";
    out << indent << indent << indent << indent << "item1[1] == item[0] or\n";
    out << indent << indent << indent << indent << "item1[0] == item[1] or\n";
    out << indent << indent << indent << indent << "item1[1] == item[1]): continue\n\n";

    out << indent << indent << "nodes = set([item1[0],item1[1],item[0],item[1]])\n";
    out << indent << indent << "if len(nodes) < 4: continue       # share same vertex\n\n";

    out << indent << indent << "(a2, b2, c2) = lineEquations[item]\n\n";

    out << indent << indent << "segment2 = [C[0] - A[0], C[1] - A[1], 0.0]\n";
    out << indent << indent << "segment3 = [D[0] - A[0], D[1] - A[1], 0.0]\n\n";

    out << indent << indent << "# check that C and D aren't on the same side of AB\n";
    out << indent << indent << "if cross(segment1, segment2)[2] * cross(segment1, segment3)[2] > 0.0: continue\n\n";

    out << indent << indent << "if cross(segment1, segment2)[2] < 0.0 or cross(segment1, segment3)[2] > 0.0:\n";
    out << indent << indent << indent << "(C, D, a2, b2, c2) = (D, C, -a2, -b2, -c2)\n";
    out << indent << indent << indent << "item = (item[1], item[0])\n\n";

    out << indent << indent << "denominator = a1 * b2 - b1 * a2\n";
    out << indent << indent << "if math.fabs(denominator) < 1e-9: continue\n\n";

    out << indent << indent << "px = (b1 * c2 - b2 * c1) / denominator\n";
    out << indent << indent << "py = (a2 * c1 - a1 * c2) / denominator\n\n";

    out << indent << indent << "distanceCD = distance(C, D)\n\n";

    out << indent << indent << "distanceAP = distance(A, [px, py])\n";
    out << indent << indent << "distanceBP = distance(B, [px, py])\n";
    out << indent << indent << "distanceCP = distance(C, [px, py])\n";
    out << indent << indent << "distanceDP = distance(D, [px, py])\n\n";

    out << indent << indent << "# check intersection is between AB and CD\n";
    out << indent << indent << "check1 = (distanceAP - 1e-6 < distanceAB and distanceBP - 1e-6 < distanceAB)\n";
    out << indent << indent << "check2 = (distanceCP - 1e-6 < distanceCD and distanceDP - 1e-6 < distanceCD)\n";
    out << indent << indent << "if not check1 or not check2: continue\n\n";

    out << indent << indent << "if math.fabs(dot(segment1, [D[0] - C[0], D[1] - C[1], 0.0]) / (distanceAB * distanceCD) + 1.0) < 1e-9: continue\n\n";

    out << indent << indent << "if ((distanceAP < 1e-6) or\n";
    out << indent << indent << indent << "(distanceBP < 1e-6 and distanceCP < 1e-6) or\n";
    out << indent << indent << indent << "(distanceDP < 1e-6)):\n";
    out << indent << indent << indent << "continue\n\n";

    out << indent << indent << "direction = [D[0] - C[0], D[1] - C[1]]\n";
    out << indent << indent << "direction = [dot(matrix[0], direction), dot(matrix[1], direction)]\n";
    out << indent << indent << "angle = math.fmod(math.atan2(direction[1], direction[0]) + 2 * math.pi, 2 * math.pi)\n\n";

    out << indent << indent << "newSegment = ([px, py], item[1], distanceAP, angle)\n\n";

    out << indent << indent << "if distanceCP < 1e-6 and angle > origMinAngle: continue\n";
    out << indent << indent << "if distanceBP > 1e-6 and angle > math.pi: continue\n\n";

    out << indent << indent << "if len(intersectingSegments) == 0:\n";
    out << indent << indent << indent << "intersectingSegments.append(newSegment)\n";
    out << indent << indent << indent << "minAngle = angle\n";
    out << indent << indent << "else:\n";
    out << indent << indent << indent << "if intersectingSegments[0][2] - 1e-9 > distanceAP:\n";
    out << indent << indent << indent << indent << "intersectingSegments[0] = newSegment\n";
    out << indent << indent << indent << indent << "minAngle = angle\n";
    out << indent << indent << indent << "elif intersectingSegments[0][2] + 1e-9 > distanceAP and angle < minAngle:\n";
    out << indent << indent << indent << indent << "intersectingSegments[0] = newSegment\n";
    out << indent << indent << indent << indent << "minAngle = angle\n\n";

    out << indent << "return intersectingSegments\n\n";

    out << "def projectToPlane(normal):\n";
    out << indent << "fh = open('" << fname << "_ElementPositions.gpl','w')\n";
    // out << indent << "fg = open('test2.dat', 'w')\n";
    out << indent << "normal = normalize(normal)\n";
    out << indent << "ori = getQuaternion(normal, [0, 0, 1])\n\n";

    out << indent << "left2Right = [0, 0, 1]\n";
    out << indent << "if math.fabs(math.fabs(dot(normal, left2Right)) - 1) < 1e-9:\n";
    out << indent << indent << "left2Right = [1, 0, 0]\n";
    out << indent << "rotL2R = ori.rotate(left2Right)\n";
    out << indent << "deviation = math.atan2(rotL2R[1], rotL2R[0])\n";
    out << indent << "rotBack = Quaternion([math.cos(0.5 * deviation), 0, 0, -math.sin(0.5 * deviation)])\n";
    out << indent << "ori = rotBack * ori\n\n";

    out << indent << "decodeVertices()\n\n";

    out << indent << "for i in range(len(vertices)):\n";
    out << indent << indent << "positions = {}\n";
    out << indent << indent << "connections = {}\n";
    out << indent << indent << "for j in range(0, len(vertices[i]), 3):\n";
    out << indent << indent << indent << "nextPos3D = ori.rotate(vertices[i][j:j+3])\n";
    out << indent << indent << indent << "nextPos2D = nextPos3D[0:2]\n";
    out << indent << indent << indent << "positions[j/3] = nextPos2D\n\n";

    out << indent << indent << "if len(positions) == 0:\n";
    out << indent << indent << indent << "continue\n\n";

    out << indent << indent << "idx = 0\n";
    out << indent << indent << "maxX = positions[0][0]\n";
    out << indent << indent << "for j in positions.keys():\n";
    out << indent << indent << indent << "if positions[j][0] > maxX:\n";
    out << indent << indent << indent << indent << "maxX = positions[j][0]\n";
    out << indent << indent << indent << indent << "idx = j\n";
    out << indent << indent << indent << "if positions[j][0] == maxX and positions[j][1] > positions[idx][1]:\n";
    out << indent << indent << indent << indent << "idx = j\n\n";

    out << indent << indent << "for j in range(0, len(triangles[i]), 3):\n";
    out << indent << indent << indent << "for k in range(0, 3):\n";
    out << indent << indent << indent << indent << "vertIdx = triangles[i][j + k]\n";
    out << indent << indent << indent << indent << "if not vertIdx in connections:\n";
    out << indent << indent << indent << indent << indent << "connections[vertIdx] = []\n";
    out << indent << indent << indent << indent << "for l in range(1, 3):\n";
    out << indent << indent << indent << indent << indent << "connections[vertIdx].append(triangles[i][j + ((k + l) % 3)])\n\n";

    out << indent << indent << "numConnections = 0\n";
    out << indent << indent << "for j in connections.keys():\n";
    out << indent << indent << indent << "connections[j] = list(set(connections[j]))\n";
    out << indent << indent << indent << "numConnections += len(connections[j])\n\n";
    out << indent << indent << "numConnections /= 2\n";

    out << indent << indent << "idChanges = squashVertices(positions, connections)\n\n";

    out << indent << indent << "lineEquations = computeLineEquations(positions, connections)\n\n";

    // out << indent << indent << "for j in positions.keys():\n";
    // out << indent << indent << indent << "for k in connections[j]:\n";
    // out << indent << indent << indent << indent << "if j < k:\n";
    // out << indent << indent << indent << indent << indent << "fg.write(\"%.8f    %.8f    \\n%.8f    %.8f\\n#    %d    %d\\n\\n\" %(positions[j][0], positions[j][1], positions[k][0], positions[k][1], j, k))\n";
    // out << indent << indent << indent << "fg.write(\"\\n\")\n\n";

    out << indent << indent << "idx = idChanges[int(idx)]\n";
    out << indent << indent << "fh.write(\"%.6f    %.6f\\n\" % (positions[idx][0], positions[idx][1]))\n\n";

    out << indent << indent << "curAngle = math.pi\n";
    out << indent << indent << "origin = idx\n";
    out << indent << indent << "count = 0\n";
    out << indent << indent << "passed = []\n";
    out << indent << indent << "while (count == 0 or distance(positions[origin], positions[idx]) > 1e-9) and not count > numConnections:\n";
    out << indent << indent << indent << "nextGen = computeMinAngle(idx, curAngle, positions, connections, False)\n";
    out << indent << indent << indent << "nextIdx = nextGen[0]\n";
    out << indent << indent << indent << "direction = [positions[nextIdx][0] - positions[idx][0],\n";
    out << indent << indent << indent << indent << "         positions[nextIdx][1] - positions[idx][1]]\n";
    out << indent << indent << indent << "curAngle = math.atan2(direction[1], direction[0])\n\n";

    out << indent << indent << indent << "intersections = checkPossibleSegmentIntersection((idx, nextIdx), positions, lineEquations, nextGen[1], [])\n";
    out << indent << indent << indent << "if len(intersections) > 0:\n";
    out << indent << indent << indent << indent << "while len(intersections) > 0 and not count > numConnections:\n";
    out << indent << indent << indent << indent << indent << "fh.write(\"%.6f    %.6f\\n\" %(intersections[0][0][0], intersections[0][0][1]))\n";
    out << indent << indent << indent << indent << indent << "count += 1\n";
    out << indent << "\n";
    out << indent << indent << indent << indent << indent << "idx = intersections[0][1]\n";
    out << indent << indent << indent << indent << indent << "direction = [positions[idx][0] - intersections[0][0][0],\n";
    out << indent << indent << indent << indent << indent << "             positions[idx][1] - intersections[0][0][1]]\n";
    out << indent << indent << indent << indent << indent << "curAngle = math.atan2(direction[1], direction[0])\n";
    out << indent << indent << indent << indent << indent << "nextGen = computeMinAngle(idx, curAngle, positions, connections, False)\n";
    out << indent << "\n";
    out << indent << indent << indent << indent << indent << "intersections = checkPossibleSegmentIntersection((None, idx), positions, lineEquations, nextGen[1], intersections[0][0])\n";
    out << indent << indent << indent << "else:\n";
    out << indent << indent << indent << indent << "idx = nextIdx\n";
    out << indent << "\n";
    out << indent << indent << indent << "fh.write(\"%.6f    %.6f\\n\" % (positions[idx][0], positions[idx][1]))\n";
    out << indent << "\n";
    out << indent << indent << indent << "if idx in passed:\n";
    out << indent << indent << indent << indent << "direction1 = [positions[idx][0] - positions[passed[-1]][0],\n";
    out << indent << indent << indent << indent << indent << "          positions[idx][1] - positions[passed[-1]][1]]\n";
    out << indent << indent << indent << indent << "direction2 = [positions[origin][0] - positions[passed[-1]][0],\n";
    out << indent << indent << indent << indent << indent << "          positions[origin][1] - positions[passed[-1]][1]]\n";
    out << indent << indent << indent << indent << "dist1 = distance(positions[idx], positions[passed[-1]])\n";
    out << indent << indent << indent << indent << "dist2 = distance(positions[origin], positions[passed[-1]])\n";
    out << indent << indent << indent << indent << "if dist1 * dist2 > 0.0 and math.fabs(math.fabs(dot(direction1, direction2) / (dist1 * dist2)) - 1.0) > 1e-9:\n";
    out << indent << indent << indent << indent << indent << "sys.stderr.write(\"error: projection cycling on element id: %d, vertex id: %d\\n\" %(i, idx))\n";
    out << indent << indent << indent << indent << "break\n\n";

    out << indent << indent << indent << "passed.append(idx)\n";
    out << indent << indent << indent << "count += 1\n";
    out << indent << indent << "fh.write(\"\\n\")\n\n";

    out << indent << indent << "if count > numConnections:\n";
    out << indent << indent << indent << "sys.stderr.write(\"error: projection cycling on element id: %d\\n\" % i)\n\n";

    out << indent << indent << "for j in range(0, len(decoration[i]), 6):\n";
    out << indent << indent << indent << "for k in range(j, j + 6, 3):\n";
    out << indent << indent << indent << indent << "nextPos3D = ori.rotate(decoration[i][k:k+3])\n";
    out << indent << indent << indent << indent << "fh.write(\"%.6f    %.6f\\n\" % (nextPos3D[0], nextPos3D[1]))\n";
    out << indent << indent << indent << "fh.write(\"\\n\")\n\n";

    out << indent << "fh.close()\n\n";

    out << "if __name__ == \"__main__\":\n";
    out << indent << "parser = argparse.ArgumentParser()\n";
    out << indent << "parser.add_argument('--export-vtk', action='store_true')\n";
    out << indent << "parser.add_argument('--export-web', action='store_true')\n";
    out << indent << "parser.add_argument('--background', nargs=3, type=float)\n";
    out << indent << "parser.add_argument('--project-to-plane', action='store_true')\n";
    out << indent << "parser.add_argument('--normal', nargs=3, type=float)\n";
    out << indent << "args = parser.parse_args()\n\n";

    out << indent << "if (args.export_vtk):\n";
    out << indent << indent << "exportVTK()\n";
    out << indent << indent << "sys.exit()\n\n";

    out << indent << "if (args.export_web):\n";
    out << indent << indent << "bgcolor = []\n";
    out << indent << indent << "if (args.background):\n";
    out << indent << indent << indent << "validBackground = True\n";
    out << indent << indent << indent << "for comp in bgcolor:\n";
    out << indent << indent << indent << indent << "if comp < 0.0 or comp > 1.0:\n";
    out << indent << indent << indent << indent << indent << "validBackground = False\n";
    out << indent << indent << indent << indent << indent << "break\n";
    out << indent << indent << indent << "if (validBackground):\n";
    out << indent << indent << indent << indent << "bgcolor = args.background\n";
    out << indent << indent << "exportWeb(bgcolor)\n";
    out << indent << indent << "sys.exit()\n\n";

    out << indent << "if (args.project_to_plane):\n";
    out << indent << indent << "normal = [0.0, 1.0, 0.0]\n";
    out << indent << indent << "if (args.normal):\n";
    out << indent << indent << indent << "normal = args.normal\n";
    out << indent << indent << "projectToPlane(normal)\n";
    out << indent << indent << "sys.exit()\n\n";

    out << indent << "parser.print_help()\n";
}

MeshData MeshGenerator::getCylinder(double length,
                                    double minor,
                                    double major,
                                    double formFactor,
                                    const unsigned int numSegments) {
    double angle = 0;
    double dAngle = Physics::two_pi / numSegments;

    MeshData mesh;
    mesh.vertices_m.push_back(Vector_t(0.0));
    for (unsigned int i = 0; i < numSegments; ++ i, angle += dAngle) {
        Vector_t node(major * cos(angle), minor * sin(angle), 0);
        mesh.vertices_m.push_back(node);

        unsigned int next = (i + 1) % numSegments;
        Vektor<unsigned int, 3> baseTriangle(0u, next + 1, i + 1);
        mesh.triangles_m.push_back(baseTriangle);

        Vektor<unsigned int, 3> sideTriangle1(i + 1, next + 1, i + numSegments + 2);
        mesh.triangles_m.push_back(sideTriangle1);

        Vektor<unsigned int, 3> sideTriangle2(next + 1, next + numSegments + 2, i + numSegments + 2);
        mesh.triangles_m.push_back(sideTriangle2);
    }

    mesh.vertices_m.push_back(Vector_t(0.0, 0.0, length));
    for (unsigned int i = 0; i < numSegments; ++ i, angle += dAngle) {
        Vector_t node(formFactor * major * cos(angle), formFactor * minor * sin(angle), length);
        mesh.vertices_m.push_back(node);

        unsigned int next = (i + 1) % numSegments;
        Vektor<unsigned int, 3> topTriangle(numSegments + 1, i + numSegments + 2, next + numSegments + 2);
        mesh.triangles_m.push_back(topTriangle);
    }

    mesh.decorations_m.push_back(std::make_pair(Vector_t(0.0), Vector_t(0, 0, length)));

    return mesh;
}

MeshData MeshGenerator::getBox(double length,
                               double width,
                               double height,
                               double formFactor) {

    MeshData mesh;
    mesh.vertices_m.push_back(Vector_t(width, height, 0.0));
    mesh.vertices_m.push_back(Vector_t(-width, height, 0.0));
    mesh.vertices_m.push_back(Vector_t(-width, -height, 0.0));
    mesh.vertices_m.push_back(Vector_t(width, -height, 0.0));

    mesh.vertices_m.push_back(Vector_t(formFactor * width, formFactor * height, length));
    mesh.vertices_m.push_back(Vector_t(-formFactor * width, formFactor * height, length));
    mesh.vertices_m.push_back(Vector_t(-formFactor * width, -formFactor * height, length));
    mesh.vertices_m.push_back(Vector_t(formFactor * width, -formFactor * height, length));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(0, 2, 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(0, 3, 2));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(0, 1, 4));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 1, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(1, 2, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(5, 2, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 3, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(6, 3, 7));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(3, 0, 7));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(7, 0, 4));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 5, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 6, 7));

    mesh.decorations_m.push_back(std::make_pair(Vector_t(0.0), Vector_t(0, 0, length)));

    return mesh;
}