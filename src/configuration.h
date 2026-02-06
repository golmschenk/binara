#pragma once

#include <filesystem>
#include "toml++/toml.hpp"

class Configuration
{
public:
    [[nodiscard]] bool prefix_session_directory_with_datetime() const;
    [[nodiscard]] int32_t get_number_of_threads() const;
    [[nodiscard]] std::filesystem::path get_session_directory_path() const;

private:
    [[nodiscard]] static int32_t initialize_number_of_threads(const toml::table& toml_configuration_table);
    [[nodiscard]] static std::filesystem::path initialize_session_directory_path(
        const toml::table& toml_configuration_table,
        int64_t tic_id,
        int32_t sector);
    [[nodiscard]] static std::filesystem::path initialize_data_directory_path(
        const toml::table& toml_configuration_table,
        int64_t tic_id,
        int32_t sector);
    static std::tuple<std::filesystem::path, std::filesystem::path, std::filesystem::path, std::filesystem::path,
                      std::filesystem::path, std::filesystem::path> copy_input_data_and_initialize_file_paths(
        const std::filesystem::path& data_directory_path,
        const std::filesystem::path&
        session_directory_path);
    Configuration(int64_t tic_id, int32_t sector);
    bool prefix_session_directory_with_datetime_;
    int32_t number_of_threads_;
    std::filesystem::path session_directory_path_;
    std::filesystem::path data_directory_path_;
    std::filesystem::path folded_observed_light_curve_path_;
    std::filesystem::path magnitudes_and_colors_path_;
    std::filesystem::path py_initialize_path_;
    std::filesystem::path states_path_;
    std::filesystem::path folded_observed_and_and_model_light_curves_path_;
    std::filesystem::path parameters_path_;

    friend void initialize_configuration(int64_t tic_id, int32_t sector);
};

void initialize_configuration(int64_t tic_id, int32_t sector);
Configuration& get_configuration();
