// Top-level "mie:" optics-table configuration parsing.
//
// Pure ConfigManager coverage: the aerosol type -> optics file map and directory
// must be read into ConfigData::mie in YAML declaration order (determinism), and
// the strict schema must keep accepting the "mie" key.  No Fortran bridge is
// involved, so this test needs no optics files on disk.
#include "catchem_config_manager.hpp"
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

    void write(const std::string& path, const std::string& body) {
        std::ofstream out(path);
        out << body;
        out.close();
    }

    bool threw_containing(const std::function<void()>& fn, const std::string& expected) {
        try {
            fn();
        } catch (const std::exception& error) {
            return std::string(error.what()).find(expected) != std::string::npos;
        }
        return false;
    }

} // namespace

int main() {
    int failures = 0;
    auto check = [&](bool condition, const std::string& label) {
        std::cout << (condition ? "  PASS: " : "  FAIL: ") << label << '\n';
        if (!condition)
            ++failures;
    };

    std::cout << "RUNNING TEST: Mie (optics table) configuration parsing\n";

    // 1. Full section parses with directory and ordered file pairs.
    {
        write("mie_full.yml",
              "simulation:\n  nx: 1\n  ny: 1\n  nz: 1\n"
              "mie:\n"
              "  directory: \"./ExtData/monochromatic/\"\n"
              "  files:\n"
              "    SS: optics_SS.v3_3.nc\n"
              "    DU: optics_DU.v15_3.nc\n"
              "    BC: optics_BC.v1_3.nc\n"
              "    NI: optics_NI.v2_5.nc\n");
        catchem::ConfigManager cfg;
        cfg.load_from_file("mie_full.yml");
        check(cfg.data.mie.directory == "./ExtData/monochromatic/", "mie.directory parsed");
        check(cfg.data.mie.files.size() == 4, "mie.files has four entries");
        bool order = cfg.data.mie.files.size() == 4 && cfg.data.mie.files[0].first == "SS" &&
                     cfg.data.mie.files[1].first == "DU" && cfg.data.mie.files[2].first == "BC" &&
                     cfg.data.mie.files[3].first == "NI";
        check(order, "mie.files preserves YAML declaration order");
        bool values = cfg.data.mie.files.size() == 4 && cfg.data.mie.files[0].second == "optics_SS.v3_3.nc" &&
                      cfg.data.mie.files[3].second == "optics_NI.v2_5.nc";
        check(values, "mie.files values are the optics file names");
        // The strict schema allow-lists "mie": no schema issue may name it.
        // (The minimal YAML legitimately fails mechanism validation, which is
        // unrelated to the mie section.)
        {
            const auto& report = cfg.validate(true);
            bool mie_rejected = false;
            for (const auto& issue : report.issues)
                if (issue.category == "schema" && issue.path == "mie")
                    mie_rejected = true;
            check(!mie_rejected, "strict schema accepts the top-level mie key");
        }
    }

    // 2. Absent section keeps the defaults (no optics configured).
    {
        write("mie_absent.yml", "simulation:\n  nx: 1\n  ny: 1\n  nz: 1\n");
        catchem::ConfigManager cfg;
        cfg.load_from_file("mie_absent.yml");
        check(cfg.data.mie.directory == "./", "absent mie section keeps default directory './'");
        check(cfg.data.mie.files.empty(), "absent mie section keeps empty file list");
    }

    // 3. Directory omitted -> default './' with files still parsed.
    {
        write("mie_nodefault.yml",
              "simulation:\n  nx: 1\n  ny: 1\n  nz: 1\n"
              "mie:\n"
              "  files:\n"
              "    SU: optics_SU.v1_3.nc\n");
        catchem::ConfigManager cfg;
        cfg.load_from_file("mie_nodefault.yml");
        check(cfg.data.mie.directory == "./", "mie.directory omitted defaults to './'");
        check(cfg.data.mie.files.size() == 1 && cfg.data.mie.files[0].first == "SU", "single file entry parsed");
    }

    // 4. Empty file value is rejected at load (fail loud, FR-009 adjacent).
    {
        write("mie_emptyfile.yml",
              "simulation:\n  nx: 1\n  ny: 1\n  nz: 1\n"
              "mie:\n"
              "  files:\n"
              "    DU: \"\"\n");
        catchem::ConfigManager cfg;
        check(threw_containing([&] { cfg.load_from_file("mie_emptyfile.yml"); }, "no optics file name"),
              "empty optics file name is rejected naming the type");
    }

    std::cout << (failures == 0 ? "SUCCESS: all mie configuration assertions passed.\n"
                                : "FAILURE: " + std::to_string(failures) + " assertion(s) failed.\n");
    return failures == 0 ? 0 : 1;
}
