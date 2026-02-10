#include "configuration.h"
#include "eclipse_method_enum.h"

#include <iostream>
#include <memory>
#include <format>
#include <omp.h>
#include <unordered_set>
#include <chrono>
#include <string>

#include "toml++/toml.hpp"

namespace
{
    std::unique_ptr<Configuration> configuration_instance = nullptr;

    const std::unordered_set<std::string> allowed_paths = {
        "output.prefix_session_directory_with_datetime",
        "system.number_of_threads",
        "modeling.use_g_magnitude",
        "modeling.use_colors",
        "modeling.use_secular_drift",
        "modeling.eclipse_method",
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
                if (!allowed_paths.contains(subpath.str()))
                {
                    std::cerr << "Unexpected configuration option in `" << subpath.str() <<
                        "` file: `" << subpath.str() << "`. Halting program.\n";
                    exit(100);
                }
            }
        }
    }

    std::string get_datetime_string()
    {
        const auto now = std::chrono::system_clock::now();
        const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        const std::tm* local_datetime = std::localtime(&now_time);

        std::ostringstream stream;
        stream << std::put_time(local_datetime, "%Y_%m_%d_%H_%M_%S");
        return stream.str();
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

bool Configuration::initialize_should_use_g_magnitude(const toml::table& toml_configuration_table)
{
    const bool should_use_g_magnitude = toml_configuration_table.at_path("modeling.use_g_magnitude").value_or(true);
    return should_use_g_magnitude;
}

bool Configuration::initialize_should_use_colors(const toml::table& toml_configuration_table)
{
    const bool should_use_colors = toml_configuration_table.at_path("modeling.use_colors").value_or(false);
    return should_use_colors;
}

bool Configuration::initialize_should_use_secular_drift(const toml::table& toml_configuration_table)
{
    const bool should_secular_drift = toml_configuration_table.at_path("modeling.use_secular_drift").value_or(false);
    return should_secular_drift;
}

EclipseMethod Configuration::initialize_eclipse_method(const toml::table& toml_configuration_table)
{
    const std::string eclipse_method_string = std::string(toml_configuration_table.at_path(
        "modeling.eclipse_method").value_or("limb_darkening"));
    EclipseMethod eclipse_method;
    if (eclipse_method_string == std::string_view("off"))
    {
        eclipse_method = EclipseMethod::OFF;
    }
    else if (eclipse_method_string == std::string_view("regular"))
    {
        eclipse_method = EclipseMethod::REGULAR;
    }
    else if (eclipse_method_string == std::string_view("limb_darkening"))
    {
        eclipse_method = EclipseMethod::LIMB_DARKENING;
    }
    else
    {
        std::cerr << "Invalid eclipse method `" << eclipse_method_string
            << "` specified in configuration.toml. Halting program.";
        exit(106);
    }
    return eclipse_method;
}

std::filesystem::path Configuration::initialize_session_directory_path(const int64_t tic_id, const int32_t sector) const
{
    std::filesystem::path sessions_root_directory{"sessions"};
    auto session_name = std::format("tic_id_{}_sector_{}", tic_id, sector);
    if (prefix_session_directory_with_datetime_)
    {
        session_name = std::format("{}_{}", get_datetime_string(), session_name);
    }
    std::filesystem::path session_directory = sessions_root_directory / session_name;
    if (std::filesystem::exists(session_directory))
    {
        std::cerr << "Output session directory `" << session_directory << "` already exists. Halting program.";
        exit(105);
    }
    std::filesystem::create_directories(session_directory);
    return session_directory;
}

std::filesystem::path Configuration::initialize_input_data_directory_path(
    const int64_t tic_id,
    const int32_t sector)
{
    std::filesystem::path input_data_directory{std::format("input_data/tic_id_{}_sector_{}", tic_id, sector)};
    if (!std::filesystem::exists(input_data_directory))
    {
        std::cerr << "Expected input data directory `" << input_data_directory << "` not found. Halting program.";
        exit(104);
    }
    return input_data_directory;
}

std::tuple<std::filesystem::path, std::filesystem::path, std::filesystem::path, std::filesystem::path,
           std::filesystem::path, std::filesystem::path> Configuration::copy_input_data_and_initialize_file_paths(
    const std::filesystem::path& data_directory_path, const std::filesystem::path& session_directory_path)
{
    if (std::filesystem::exists("configuration.toml"))
    {
        std::filesystem::copy("configuration.toml", session_directory_path / "configuration.toml");
    }
    std::filesystem::path folded_observed_light_curve_path = session_directory_path / "folded_observed_light_curve.txt";
    std::filesystem::copy(data_directory_path / "folded_observed_light_curve.txt", folded_observed_light_curve_path);
    std::filesystem::path magnitudes_and_colors_path = session_directory_path / "magnitudes_and_colors.txt";
    std::filesystem::copy(data_directory_path / "magnitudes_and_colors.txt", magnitudes_and_colors_path);
    std::filesystem::path py_initialize_path = session_directory_path / "py_initialize.txt";
    std::filesystem::copy(data_directory_path / "py_initialize.txt", py_initialize_path);
    std::filesystem::path states_path = session_directory_path / "states.txt";
    std::ofstream states_ofstream(states_path);
    std::filesystem::path folded_observed_and_and_model_light_curves_path = session_directory_path /
        "folded_observed_and_model_light_curves.txt";
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
    session_directory_path_ = initialize_session_directory_path(tic_id, sector);
    input_data_directory_path_ = initialize_input_data_directory_path(tic_id, sector);
    std::tie(folded_observed_light_curve_path_, magnitudes_and_colors_path_, py_initialize_path_, states_path_,
             folded_observed_and_and_model_light_curves_path_,
             parameters_path_) = copy_input_data_and_initialize_file_paths(
        input_data_directory_path_, session_directory_path_);
    should_use_g_magnitude_ = initialize_should_use_g_magnitude(toml_configuration_table);
    should_use_colors_ = initialize_should_use_colors(toml_configuration_table);
    should_use_secular_drift_ = initialize_should_use_secular_drift(toml_configuration_table);
    eclipsing_method_ = initialize_eclipse_method(toml_configuration_table);
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

std::filesystem::path Configuration::get_py_initialize_path() const
{
    return py_initialize_path_;
}

std::filesystem::path Configuration::get_states_path() const
{
    return states_path_;
}

std::filesystem::path Configuration::get_folded_observed_and_model_light_curves_path() const
{
    return folded_observed_and_and_model_light_curves_path_;
}

std::filesystem::path Configuration::get_parameters_path() const
{
    return parameters_path_;
}

bool Configuration::should_use_g_magnitude() const
{
    return should_use_g_magnitude_;
}

bool Configuration::should_use_colors() const
{
    return should_use_colors_;
}

bool Configuration::should_use_secular_drift() const
{
    return should_use_secular_drift_;
}

EclipseMethod Configuration::eclipsing_method() const
{
    return eclipsing_method_;
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
