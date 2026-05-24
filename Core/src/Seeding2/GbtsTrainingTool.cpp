#include<iostream>
#include <string>
#include <fstream>

#include "Acts/Seeding2/GbtsTrainingTool.hpp"


namespace Acts::Experimental {

LayerDescription::LayerDescription(float minR_, float maxR_, float minZ_, float maxZ_, std::int32_t gbtsId_) :
    minR(minR_), maxR(maxR_), minZ(minZ_), maxZ(maxZ_), gbtsId(gbtsId_){}


GbtsTrainingTool::GbtsTrainingTool(std::istream& inStream){

    //define how many lines there are
    std::uint32_t lines{};
    std::string line{};
    while(std::getline(inStream, line)){
        lines++;
    }
    inStream.clear();
    inStream.seekg(0);

    m_layerGeometry.reserve(lines);

    // create geometry objects
    float minR{};
    float maxR{};
    

    float minZ{};
    float maxZ{};

    std::int32_t gbtsId{};

    for(std::uint32_t l = 0; l < lines; l++) {

        inStream >> minR >> maxR
                 >> minZ >> maxZ
                 >> gbtsId;

        m_layerGeometry.emplace_back(minR, maxR, minZ, maxZ, gbtsId);
    }

    // create linked pairs using vector of indices to m_layerGeometry
    for(std::uint32_t i = 0; i < m_layerGeometry.size(); i++){
        for(std::uint32_t j = i + 1; j < m_layerGeometry.size(); j++){

            LayerIdPair pair;
            pair.first = m_layerGeometry[i].gbtsId;
            pair.second = m_layerGeometry[j].gbtsId;

            // key = GBTS id, value = number of transitions 
            m_layerPairs.emplace(pair, 0);
        }
    }

} 

void GbtsTrainingTool::addTrack(const std::vector<TrackCoordinates>& track){

    if(track.size() < 2){

        std::cout<< "Warning: Track only has one measurement"<<std::endl;
        return;
    }
    //container for gbts IDs of the track
    std::vector<std::int32_t> layerGbtsIds{};
    layerGbtsIds.reserve(track.size());
    
    // find GBTS ids for all measurements in a track
    for(const auto& measurement : track){

        const float r = measurement.r;
        const float z = measurement.z;

        const std::int32_t gbtsId = findGbtsIdByCoord(r, z);

        layerGbtsIds.emplace_back(gbtsId);
    }

    // update map with track layer transitions
    for(std::uint32_t id = 0; id + 1 < layerGbtsIds.size(); id++){

        const std::int32_t index1 = layerGbtsIds[id];
        const std::int32_t index2 = layerGbtsIds[id + 1];

        if(index1 == index2){
            std::cout << "Warning: track transitions between same layer"<<std::endl;

            continue;
        }

        m_layerPairs[{index1, index2}] += 1;
    }

    m_totalTracks += 1;

}   

void GbtsTrainingTool::createConnectionTable(const std::filesystem::path& outputFileLocation, const double probThreshold) const {

    if (m_totalTracks == 0) {
    std::runtime_error("Warning: no tracks were added when creating connection table");
    return;
}
    //define output text file
    std::ofstream outputFile(outputFileLocation);

    // for every layer transition, only output those that pass probability cut
    for(const auto& [indexes, nTransitions] : m_layerPairs){

        const double probability = static_cast<double>(nTransitions) / static_cast<double>(m_totalTracks);

        if(probThreshold < probability){
            
            outputFile << indexes.first << " " << indexes.second << "\n";
        }
    }
}

std::int32_t GbtsTrainingTool::findGbtsIdByCoord(float r, float z) const {

    for(const auto& layer : m_layerGeometry){
        
        if(layer.minZ < z && z < layer.maxZ){
            if(layer.minR < r && r < layer.maxR){

                return layer.gbtsId;
            }
        }
    }

    throw std::runtime_error("No matching GBTS ID found for coordinates");
}

} //namespace Acts::Experimental