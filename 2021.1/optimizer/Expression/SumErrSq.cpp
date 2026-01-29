#include "Expression/SumErrSq.h"

const std::string SumErrSq::name = "SumErrSq";

std::istream & operator>>(std::istream & stream, Measurement & measurement) {
    stream >> measurement.spos >> measurement.measurement;
    return stream;
}

/**
 *  parses a simple list of spos and measurements using tab as delimiter.
 */
void SumErrSq::parseMeasurements(std::string measurement_filename) {

    measurements_.clear();

    std::ifstream measurements_file;
    measurements_file.open(measurement_filename.c_str(), std::ios::in);
    if(!measurements_file) {
        throw OptPilotException("SumErrSq::parseMeasurements()",
                "Error opening file " + measurement_filename);
    }

    std::copy(std::istream_iterator<Measurement>(measurements_file),
              std::istream_iterator<Measurement>(),
              std::back_inserter(measurements_));

    measurements_file.close();
}

