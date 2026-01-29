namespace interpolation {

std::ostream& operator << (std::ostream& out,
                           const PolynomialCoefficient& coeff) {
    std::vector<int> in_var = coeff.InVariables();
    for (size_t i = 0; i < coeff.InVariables().size(); ++i)
        out << in_var[i] << " ";
    out << "| " << coeff.OutVariable();
    out << " | " << coeff.Coefficient();
    return out;
}

void PolynomialCoefficient::SpaceTransform(std::vector<int> space_in,
                                           std::vector<int> space_out) {
    std::map<int, int> mapping; //probably optimise this
    for (unsigned int i = 0; i < space_out.size(); i++)
        for (unsigned int j = 0; j < space_in.size(); j++)
            if (space_out[i] == space_in[j])
                mapping[j] = i; //mapping[space_in_index] returns space_out_index

    std::vector<int> in_variables(_inVarByVec.size());
    for (unsigned int con=0; con<in_variables.size(); con++)
    {
        if (mapping.find(in_variables[con]) != mapping.end())
           in_variables[con] = mapping[in_variables[con]];
      else
          throw(GeneralClassicException(
                "PolynomialVector::PolynomialCoefficient::SpaceTransform",
                "Input variable not found in space transform"
          ));
    }

    if(mapping.find(_outVar) != mapping.end())
        _outVar = mapping[_outVar];
    else
        throw(GeneralClassicException(
              "PolynomialVector::PolynomialCoefficient::SpaceTransform",
              "Output variable not found in space transform"
        ));
    _inVarByVec = in_variables;
}
}

