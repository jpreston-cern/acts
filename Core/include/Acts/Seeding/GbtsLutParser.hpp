#pragma once

#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <array>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>





namespace Acts::Experimental{

class GbtsLutParser{

    public:
        explicit GbtsLutParser(SeedFinderGbtsConfig& config){

            if (config.useML) {
    if (config.lutInputFile.empty()) {
      throw std::runtime_error("Cannot find ML predictor LUT file");
    } else {
      m_mlLUT.reserve(100);
      std::ifstream ifs(std::string(config.lutInputFile).c_str());

      if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open LUT file");
      }

      float cl_width{}, min1{}, max1{}, min2{}, max2{};

      while (ifs >> cl_width >> min1 >> max1 >> min2 >> max2) {
        std::array<float, 5> lut_line = {cl_width, min1, max1, min2, max2};
        m_mlLUT.emplace_back(lut_line);
      }
      if (!ifs.eof()) {
        // ended if parse error present, not clean EOF

        throw std::runtime_error("Stopped reading LUT file due to parse error");
      }

      ifs.close();
    }
  }
}

std::vector<std::array<float, 5>>& mlLUT(){
    
    return m_mlLUT;
}

private:
    std::vector<std::array<float, 5>> m_mlLUT{};

};

}//name space Acts::Exerpimental