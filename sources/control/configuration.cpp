//
// Created by 159-egi on 2/16/17.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>

Configuration::Configuration(json11::Json config_part) :
    m_config(std::move(config_part)) {
}

Configuration::Configuration(string filename) {
    if (filename.empty()) {
        throw runtime_error("Configuration filename is not given.");
    }

    for(int r = 0; r < mpi::size(); ++r) {
        mpi::barrier();
        if (r != mpi::rank()) {
            continue;
        }
        std::ifstream config_dc(filename, std::ios::in);

        if (config_dc) {
            string errors;
            string config_str = string(std::istreambuf_iterator<char>(config_dc),
                                       std::istreambuf_iterator<char>());

            m_config = json11::Json::parse(config_str, errors);

            if (!errors.empty()) {
                throw runtime_error("ERROR! Config syntax is incorrect: " + errors);
            }
        } else {
            throw runtime_error("ERROR! Cant open configuration file: " + filename);
        }

        config_dc.close();
    }
}

bool Configuration::exist(const string& name) const {
    return !m_config[name].is_null();
}

bool Configuration::exist(const string& section,
                          const string& name) const {
    return !m_config[section][name].is_null();
}

bool Configuration::exist(const string& section,
                          const string& subsection,
                          const string& name) const {
    return !m_config[section][subsection][name].is_null();
}

Configuration Configuration::operator()(const string &name) const {
    if (m_config[name].is_null()) {
        throw runtime_error("Can't find parameter '" + name + "' in configuration.");
    }
    return Configuration(m_config[name]);
}

Configuration Configuration::operator()(const string &section,
                                        const string &name) const {
    if (m_config[section].is_null() ||
        m_config[section][name].is_null()) {

        throw runtime_error("Can't find parameter '" + name +
                            "' in section '" + section + "'.");
    }
    return Configuration(m_config[section][name]);
}

Configuration Configuration::operator()(const string &section,
                                        const string &subsection,
                                        const string &name) const {
    if (m_config[section].is_null() ||
        m_config[section][subsection].is_null() ||
        m_config[section][subsection][name].is_null()) {

        throw runtime_error("Can't find parameter '" + name +
                                    "' in subsection '" + subsection +
                                    "' of section '" + section + "'.");
    }
    return Configuration(m_config[section][subsection][name]);
}

bool Configuration::is_array() const {
    return m_config.is_array();
}

bool Configuration::to_bool() const {
    if (!m_config.is_bool()) {
        throw runtime_error("Can't parse boolean from configuration.");
    }
    return m_config.bool_value();
}

int Configuration::to_int() const {
    if (!m_config.is_number()) {
        throw runtime_error("Can't parse int from configuration.");
    }
    return m_config.int_value();
}

uint Configuration::to_uint() const {
    if (!m_config.is_number() || m_config.int_value() < 0) {
        throw runtime_error("Can't parse uint from configuration.");
    }
    return static_cast<uint>(m_config.int_value());
}

double Configuration::to_double() const {
    if (!m_config.is_number()) {
        throw runtime_error("Can't parse double from configuration.");
    }
    return m_config.number_value();
}

string Configuration::to_string() const {
    if (!m_config.is_string()) {
        throw runtime_error("Can't parse string from configuration.");
    }
    return m_config.string_value();
}

vector<Configuration> Configuration::to_array() const {
    vector<Configuration> elements;
    for (auto &element: m_config.array_items()) {
        elements.emplace_back(Configuration(element));
    }
    return elements;
}

std::map<string, double> Configuration::to_map() const {
    std::map<string, double> elements;
    if (!m_config.is_object()) {
        throw runtime_error("Configuration::to_map() error: Not an object");
    }

    for (auto& element: m_config.object_items()) {
        if (element.second.is_number()) {
            elements[element.first] = element.second.number_value();
        }
        else {
            throw runtime_error("Configuration::to_map() error: Not a number");
        }
    }
    return elements;
}