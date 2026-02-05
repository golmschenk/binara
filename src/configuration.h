#pragma once

#include <filesystem>
#include "toml++/toml.hpp"

class Configuration
{
public:
    [[nodiscard]] bool prefix_session_directory_with_datetime() const;
    [[nodiscard]] int32_t number_of_threads() const;
    [[nodiscard]] std::filesystem::path session_directory() const;

private:
    [[nodiscard]] static int32_t initialize_number_of_threads(const toml::table& toml_configuration_table);
    [[nodiscard]] static std::filesystem::path initialize_session_directory(
        const toml::table& toml_configuration_table);
    Configuration();
    bool prefix_session_directory_with_datetime_;
    int32_t number_of_threads_;

    friend void initialize_configuration();
};

void initialize_configuration();
Configuration& get_configuration();
