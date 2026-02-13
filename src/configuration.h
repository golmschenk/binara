#pragma once

#include <filesystem>

#include "eclipse_method_enum.h"
#include "toml++/toml.hpp"

class Configuration
{
public:
    [[nodiscard]] bool prefix_session_directory_with_datetime() const;
    [[nodiscard]] int32_t get_number_of_threads() const;
    [[nodiscard]] std::filesystem::path get_session_directory_path() const;
    [[nodiscard]] std::filesystem::path get_py_initialize_path() const;
    [[nodiscard]] std::filesystem::path get_states_path() const;
    [[nodiscard]] std::filesystem::path get_folded_observed_and_model_light_curves_path() const;
    [[nodiscard]] std::filesystem::path get_parameters_path() const;
    [[nodiscard]] bool should_use_g_magnitude() const;
    [[nodiscard]] bool should_use_colors() const;
    [[nodiscard]] bool should_use_secular_drift() const;
    [[nodiscard]] EclipseMethod eclipsing_method() const;

private:
    [[nodiscard]] static int32_t initialize_number_of_threads(const toml::table& toml_configuration_table);
    [[nodiscard]] static bool initialize_should_use_g_magnitude(const toml::table& toml_configuration_table);
    [[nodiscard]] static bool initialize_should_use_colors(const toml::table& toml_configuration_table);
    [[nodiscard]] static bool initialize_should_use_secular_drift(const toml::table& toml_configuration_table);
    [[nodiscard]] static EclipseMethod initialize_eclipse_method(const toml::table& toml_configuration_table);
    [[nodiscard]] std::filesystem::path initialize_session_directory_path(int64_t tic_id, int32_t sector) const;
    [[nodiscard]] static std::filesystem::path initialize_input_data_directory_path(int64_t tic_id, int32_t sector);
    static std::tuple<std::filesystem::path, std::filesystem::path, std::filesystem::path, std::filesystem::path,
                      std::filesystem::path, std::filesystem::path> copy_input_data_and_initialize_file_paths(
        const std::filesystem::path& data_directory_path,
        const std::filesystem::path&
        session_directory_path);
    Configuration(int64_t tic_id, int32_t sector);
    bool prefix_session_directory_with_datetime_;
    bool should_use_g_magnitude_;
    bool should_use_colors_;
    bool should_use_secular_drift_;
    EclipseMethod eclipsing_method_;
    int32_t number_of_threads_;
    std::filesystem::path session_directory_path_;
    std::filesystem::path input_data_directory_path_;
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
