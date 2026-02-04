#include "configuration.h"

#include <iostream>
#include <unordered_set>

#include "toml++/toml.hpp"

namespace
{
    Configuration* configuration_instance = nullptr;

    const std::unordered_set<std::string> allowed_paths = {
        "output.prefix_session_directory_with_datetime",
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

Configuration::Configuration()
{
    toml::table toml_configuration_table;
    if (std::filesystem::exists("configuration.toml"))
    {
        toml_configuration_table = toml::parse_file("configuration.toml");
    }

    // Check for invalid configuration keys (to prevent typos getting through).
    validate_table(toml_configuration_table, toml::path(""));

    // Set configuration fields.
    prefix_session_directory_with_datetime_ = toml_configuration_table.at_path(
        "output.prefix_session_directory_with_datetime").value_or(true);
    std::cout << "prefix_session_directory_with_datetime = " << std::boolalpha <<
        prefix_session_directory_with_datetime_ << std::endl;
    number_of_threads_ = toml_configuration_table.at_path("system.number_of_threads").value_or(-1);
    std::cout << "prefix_session_directory_with_datetime = " << number_of_threads_ << std::endl;
}

bool Configuration::prefix_session_directory_with_datetime() const
{
    return prefix_session_directory_with_datetime_;
}

int64_t Configuration::number_of_threads() const
{
    return number_of_threads_;
}

void initialize_configuration()
{
    if (configuration_instance != nullptr)
    {
        std::cerr << "Configuration already initialized. Halting program.\n";
        exit(101);
    }

    configuration_instance = new Configuration();
}

Configuration& get_configuration()
{
    if (configuration_instance == nullptr)
    {
        std::cerr << "Configuration not initialized. Call initialize_configuration first. Halting program.\n";
        exit(102);
    }

    return *configuration_instance;
}
