//
// Class FromFile
//   This class parses a file that contains design variable values.
//   Each column belongs to a design variable.
//   The first line is considered as header and consists of the
//   design variable name. The name has to agree with the string
//   in the input file.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef SAMPLE_FROMFILE_H
#define SAMPLE_FROMFILE_H

#include "Sample/SamplingMethod.h"

#include <string>
#include <vector>

/**
 * @class FromFile
 * @brief Sampling method that reads design variable values from a text file.
 *
 * This class parses a file containing design variable samples. Each column corresponds to a variable,
 * and the first line must contain the names of the variables. It supports selecting a specific column
 * and accessing values sequentially or cyclically.
 *
 * The expected file format is:
 * @code
 * a b c
 * 1 2 3
 * 4 5 6
 * ...
 * @endcode
 * where "a", "b", "c" are variable names and subsequent lines are numerical values.
 */
class FromFile : public SamplingMethod {

public:
    explicit FromFile(const std::string& filename, const std::string& dvarName, std::size_t modulo);

    /// Destructor
    ~FromFile() override = default;

    // Disable copying to avoid issues with memory management
    FromFile(const FromFile&) = delete;
    FromFile& operator=(const FromFile&) = delete;

    // Allow move semantics for efficient resource transfer
    FromFile(FromFile&&) noexcept = default;
    FromFile& operator=(FromFile&&) noexcept = default;

    /**
     * @brief Assign a sampled value to an individual's gene.
     *
     * @param ind The individual to modify.
     * @param i Index of the gene (design variable) to assign.
     */
    void create(std::shared_ptr<SampleIndividual>& ind, std::size_t i) override;

    /**
     * @brief Parses and loads the data from the file into memory.
     *
     * @param args Command-line arguments (unused here).
     * @param comm Communication context (unused here).
     *
     * @throws OpalException if the file is invalid or the variable is missing.
     */
    void allocate(const CmdArguments_t& args, const Comm::Bundle_t& comm) override;

    /**
     * @brief Returns the next value for the given individual ID.
     *
     * The value is selected cyclically from the chain.
     *
     * @param id The individual's ID.
     * @return The corresponding sample value.
     */
    double getNext(unsigned int id);

    /**
     * @brief Get the number of lines in the file (including the header).
     *
     * @return Number of total lines parsed from the file.
     */
    unsigned int getSize() const;

private:
    /// The values for the selected design variable loaded from the file
    std::vector<double> chain_m;

    /// Modulo used to wrap indices
    std::size_t mod_m;

    /// File name where samples are read from
    std::string filename_m;

    /// Name of the design variable to extract
    std::string dvarName_m;

    /// Number of lines in the file (including header)
    std::size_t globalSize_m;
};

#endif
