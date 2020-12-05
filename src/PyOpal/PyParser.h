namespace PyOpal {
namespace parser {

std::string parser_docstring = "parser module parses the input";

std::string initialise_from_opal_file_docstring = 
  "initialise_from_opal_file(file_name)\n"
  "\n"
  "Initialise a lattice using a MAD-style opal file\n"
  "\n"
  "Parameters\n"
  "----------\n"
  "file_name : string\n"
  "    The name of the MAD-style OPAL file\n"
  "\n"
  "Returns\n"
  "-------\n"
  "Returns None\n"
  "\n"
  "Notes\n"
  "-----\n"
  "This function sometimes makes an error message from openMPI; if this is a\n"
  "problem, try rebuilding the MPI library with --disable-dlopen switch. Note\n"
  "also that, at time of writing, enabling OPAL HD5 files will crash the\n"
  "parser. Under investigation.\n";
py::object initialise_from_opal_file(std::string file_name);


std::string list_objects_docstring = 
  "list_objects(regular_expression)\n"
  "\n"
  "Make a string that lists objects that have been instantiated by OPAL\n"
  "\n"
  "Parameters\n"
  "----------\n"
  "regular_expression : string\n"
  "    A regular_expression used for matching to object names\n"
  "\n"
  "Returns\n"
  "-------\n"
  "Returns a string containing the list of objects that have been\n"
  "instantiated.\n";
std::string list_objects(std::string regular_expression);


std::string is_initialised_docstring = 
  "is_initialised()\n"
  "\n"
  "Return true if OPAL has been initialised already. OPAL can only be\n"
  "initialised once per python session.\n";
bool is_initialised();
}
}
