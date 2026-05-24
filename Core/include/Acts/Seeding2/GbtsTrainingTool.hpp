//includes needed
#include <cstdint>
#include <filesystem>
#include <istream>
#include <map>
#include <vector>

namespace Acts::Experimental {
    
struct TrackCoordinates{

    float x{};
    float y{};
    float z{};
    float r{};

};


struct LayerDescription{

    LayerDescription(float minR_, float maxR_, float minZ_, float maxZ_, std::int32_t gbtsId_);

    // r values
    float minR{};
    float maxR{};

    // z values
    float minZ{};
    float maxZ{};

    std::int32_t gbtsId{};

    
    
};

using LayerIdPair = std::pair<std::int32_t, std::int32_t>;

class GbtsTrainingTool{

    public:

    explicit GbtsTrainingTool(std::istream& inStream);

    void addTrack(const std::vector<TrackCoordinates>& Track);

    void createConnectionTable(const std::filesystem::path& outputFileLocations, const double probThreshold) const;

    private:

    std::int32_t findGbtsIdByCoord(float r, float z) const;

    std::vector<LayerDescription> m_layerGeometry{};

    std::map<LayerIdPair, std::uint32_t> m_layerPairs{};

    std::uint32_t m_totalTracks = 0;
};

}