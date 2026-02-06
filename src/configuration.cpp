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
    std::filesystem::path sessions_root_directory{"sessions"};
    std::filesystem::path session_directory = sessions_root_directory / std::format(
        "sessions/tic_id_{}_sector_{}", tic_id, sector);
    if (std::filesystem::exists(session_directory))
    {
        std::cerr << "Output session directory `" << session_directory << "` already exists. Halting program.";
        exit(105);
    }
    return session_directory;
}

std::filesystem::path Configuration::initialize_data_directory_path(
    const toml::table& toml_configuration_table,
    const int64_t tic_id,
    const int32_t sector)
{
    std::filesystem::path data_directory{std::format("data/tic_id_{}_sector_{}", tic_id, sector)};
    if (!std::filesystem::exists(data_directory))
    {
        std::cerr << "Expected data directory `" << data_directory << "` not found. Halting program.";
        exit(104);
    }
    return data_directory;
}

std::tuple<std::filesystem::path, std::filesystem::path, std::filesystem::path, std::filesystem::path,
           std::filesystem::path, std::filesystem::path> Configuration::copy_input_data_and_initialize_file_paths(
    const std::filesystem::path& data_directory_path, const std::filesystem::path& session_directory_path)
{
    if (std::filesystem::exists("configuration.toml"))
    {
        std::filesystem::copy("configuration.toml", session_directory_path / "configuration.toml");
    }
    std::filesystem::path folded_observed_light_curve_path = session_directory_path / "folded_model_light_curve.txt";
    std::filesystem::copy(data_directory_path / "folded_model_light_curve.txt", folded_observed_light_curve_path);
    std::filesystem::path magnitudes_and_colors_path = session_directory_path / "magnitudes_and_colors.txt";
    std::filesystem::copy(data_directory_path / "magnitudes_and_colors.txt", magnitudes_and_colors_path);
    std::filesystem::path py_initialize_path = session_directory_path / "py_initialize.txt";
    std::filesystem::copy(data_directory_path / "py_initialize.txt", py_initialize_path);
    std::filesystem::path states_path = session_directory_path / "states.txt";
    std::ofstream states_ofstream(states_path);
    std::filesystem::path folded_observed_and_and_model_light_curves_path = session_directory_path /
        "folded_observed_and_and_model_light_curves.txt";
    std::ofstream folded_observed_and_and_model_light_curves_ofstream(folded_observed_and_and_model_light_curves_path);
    std::filesystem::path parameters_path = session_directory_path / "parameters.txt";
    std::ofstream parameters_ofstream(parameters_path);
    return {
        folded_observed_light_curve_path, magnitudes_and_colors_path, py_initialize_path, states_path,
        folded_observed_and_and_model_light_curves_path, parameters_path
    };
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
    std::tie(folded_observed_light_curve_path_, magnitudes_and_colors_path_, py_initialize_path_, states_path_,
             folded_observed_and_and_model_light_curves_path_,
             parameters_path_) = copy_input_data_and_initialize_file_paths(
        data_directory_path_, session_directory_path_);
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
