
#pragma once

#ifdef __cplusplus
class Configuration
{
private:
    bool prefix_session_directory_with_datetime_;

    Configuration();

public:
    [[nodiscard]] bool prefix_session_directory_with_datetime() const;

    friend Configuration get_configuration();
};

Configuration get_configuration();
#endif