#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace catchem {

    struct EmissionField {
        std::string field_name;
        std::string long_name;
        std::string units;
        int nx = 0;
        int ny = 0;
        int nz = 1;
        double factors = 1.0;
        bool is_loaded = false;
        bool is_valid = false;
        bool diagnostic = false;
        bool time_interpolate = true;
        std::string interpolation_method = "bilinear";

        // Point source metadata
        int npts = 0;
        std::vector<double> lat;
        std::vector<double> lon;
        std::vector<double> stkdm;
        std::vector<double> stkht;
        std::vector<double> stktk;
        std::vector<double> stkve;
        std::vector<int> ip;
        std::vector<int> jp;
        std::vector<double> pemis;
        std::vector<double> pbot;
        std::vector<double> ptop;

        // Gridded data buffer (nx * ny * nz * n_times)
        std::vector<double> emission_data;
    };

    struct EmissionCategory {
        std::string category_name;
        std::string description;
        bool is_active = true;
        bool gridded = true;
        bool is_2d = true;
        bool diagnostic = true;
        double global_scale = 1.0;
        double topfraction = -1.0;
        std::string source_file;
        std::string format;
        std::string frequency;
        std::string regrid_method = "none";
        std::string time_interpolation = "none";
        std::string vertical_dist = "none";

        std::map<std::string, EmissionField> fields;
    };

    class EmissionManager {
    public:
        std::map<std::string, EmissionCategory> categories;
        bool is_loaded = false;

        EmissionManager() = default;

        void add_category(const EmissionCategory& category);
        bool has_category(const std::string& name) const;
        const EmissionCategory* get_category(const std::string& name) const;
    };

} // namespace catchem
