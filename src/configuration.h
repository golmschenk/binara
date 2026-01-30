#pragma once

#include <filesystem>

class Configuration
{
public:
    [[nodiscard]] bool prefix_session_directory_with_datetime() const;
    [[nodiscard]] std::filesystem::path session_directory() const;

private:
    Configuration();
    bool prefix_session_directory_with_datetime_;

    friend void initialize_configuration();
};

void initialize_configuration();
Configuration& get_configuration();
