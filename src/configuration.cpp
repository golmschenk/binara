#include "configuration.h"

#include <iostream>
#include <memory>
#include <format>
#include <omp.h>
#include <unordered_set>

#include "toml++/toml.hpp"

namespace
{
    std::unique_ptr<Configuration> configuration_instance = nullptr;

    const std::unordered_set<std::string> allowed_paths = {
        "output.prefix_session_directory_with_datetime",
        "system.number_of_threads",
    };

    void validate_table(const toml::table& table, const toml::path& path)
    {
        for (const auto& [key, node] : table)
        {
            toml::path subpath;
            if (!path.str().empty())
            {
                subpath = toml::path(std::string(path.str()) + "." + std::string(key.str()));
            }
            else
            {
                subpath = toml::path(std::string(key.str()));
            }

            if (node.is_table())
            {
                validate_table(*node.as_table(), subpath);
            }
            else
            {
                if (!allowed_paths.count(subpath.str()))
                {
                    std::cerr << "Unexpected configuration option in `" << subpath.str() <<
                        "` file: `" << subpath.str() << "`. Halting program.\n";
                    exit(100);
                }
            }
        }
    }
}

int32_t Configuration::initialize_number_of_threads(const toml::table& toml_configuration_table)
{
    int32_t number_of_threads = toml_configuration_table.at_path("system.number_of_threads").value_or(-1);
    if (number_of_threads == -1)
    {
        number_of_threads = omp_get_num_procs();
        std::cout << "number_of_threads = " << number_of_threads << " (was -1 in configuration file)" << std::endl;
    }
    else if (number_of_threads <= 0)
    {
        std::cerr << "Configuration option `number_of_threads` set to `" << number_of_threads <<
            "` but must be positive or -1. Halting program." << std::endl;
        exit(103);
    }
    else
    {
        std::cout << "number_of_threads = " << number_of_threads << std::endl;
    }
    return number_of_threads;
}

std::filesystem::path Configuration::initialize_session_directory_path(
    const toml::table& toml_configuration_table,
    const int64_t tic_id,
    const int32_t sector)
{
    std::filesystem::path session_directory{std::format("session/tic_id_{}_sector_{}", tic_id, sector)};
    return session_directory;
}

std::filesystem::path Configuration::initialize_data_directory_path(
    const toml::table& toml_configuration_table,
    const int64_t tic_id,
    const int32_t sector)
{
    std::filesystem::path data_directory{std::format("data/tic_id_{}_sector_{}", tic_id, sector)};
    return data_directory;
}

Configuration::Configuration(const int64_t tic_id, const int32_t sector)
{
    toml::table toml_configuration_table;
    if (std::filesystem::exists("configuration.toml"))
    {
        toml_configuration_table = toml::parse_file("configuration.toml");
    }

    // Check for invalid configuration keys (to prevent typos getting through).
    validate_table(toml_configuration_table, toml::path(""));

    prefix_session_directory_with_datetime_ = toml_configuration_table.at_path(
        "output.prefix_session_directory_with_datetime").value_or(true);
    std::cout << "prefix_session_directory_with_datetime = " << std::boolalpha <<
        prefix_session_directory_with_datetime_ << std::endl;

    number_of_threads_ = initialize_number_of_threads(toml_configuration_table);
    session_directory_path_ = initialize_session_directory_path(toml_configuration_table, tic_id, sector);
    data_directory_path_ = initialize_data_directory_path(toml_configuration_table, tic_id, sector);
}

bool Configuration::prefix_session_directory_with_datetime() const
{
    return prefix_session_directory_with_datetime_;
}

int32_t Configuration::get_number_of_threads() const
{
    return number_of_threads_;
}

std::filesystem::path Configuration::get_session_directory_path() const
{
    return session_directory_path_;
}

void initialize_configuration(const int64_t tic_id, const int32_t sector)
{
    configuration_instance.reset(new Configuration(tic_id, sector));
}

Configuration& get_configuration()
{
    if (configuration_instance == nullptr)
    {
        std::cerr << "Configuration not initialized. Call `initialize_configuration` first. Halting program.\n";
        exit(102);
    }

    return *configuration_instance;
}
