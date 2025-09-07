#include <cstring>
#include <zephyr/io/vtk_type.h>

namespace zephyr::io {

VtkType::VtkType() {
    m_value = VtkType::Undefined;
}

VtkType::VtkType(const ValueType& value) {
    if (value <= VtkType::Float64) {
        m_value = value;
    }
    else {
        m_value = VtkType::Undefined;
    }
}

VtkType::VtkType(const char* value) {
    if (!std::strcmp(value, "Undefined")) {
        m_value = Undefined;
    } else if (!std::strcmp(value, "Int8")) {
        m_value = Int8;
    } else if (!std::strcmp(value, "Int16")) {
        m_value = Int16;
    } else if (!std::strcmp(value, "Int32")) {
        m_value = Int32;
    } else if (!std::strcmp(value, "Int64")) {
        m_value = Int64;
    } else if (!std::strcmp(value, "UInt8")) {
        m_value = UInt8;
    } else if (!std::strcmp(value, "UInt16")) {
        m_value = UInt16;
    } else if (!std::strcmp(value, "UInt32")) {
        m_value = UInt32;
    } else if (!std::strcmp(value, "UInt64")) {
        m_value = UInt64;
    } else if (!std::strcmp(value, "Float32")) {
        m_value = Float32;
    } else if (!std::strcmp(value, "Float64")) {
        m_value = Float64;
    } else {
        m_value = Undefined;
    }
}

VtkType::VtkType(const std::string& value)
    : VtkType(value.c_str()) { }

VtkType::ValueType VtkType::value_type_conversion() const {
    return m_value;
}

const char* VtkType::const_char_conversion() const {
    switch (m_value) {
        case Int8:
            return "Int8";
        case Int16:
            return "Int16";
        case Int32:
            return "Int32";
        case Int64:
            return "Int64";
        case UInt8:
            return "UInt8";
        case UInt16:
            return "UInt16";
        case UInt32:
            return "UInt32";
        case UInt64:
            return "UInt64";
        case Float32:
            return "Float32";
        case Float64:
            return "Float64";
        default:
            return "Undefined";
    }
}

std::string VtkType::string_conversion() const {
    return const_char_conversion();
}

bool VtkType::operator==(const VtkType& value) const {
    return m_value == value;
}

bool VtkType::operator!=(const VtkType& value) const {
    return m_value != value;
}

bool VtkType::is_undefined() const {
    return m_value == VtkType::Undefined;
}

bool VtkType::is_signed() const {
    return VtkType::Int8 <= m_value && m_value <= VtkType::Int64;
}

bool VtkType::is_unsigned() const {
    return VtkType::UInt8 <= m_value && m_value <= VtkType::UInt64;
}

bool VtkType::is_floating() const {
    return m_value == VtkType::Float32 || m_value == VtkType::Float64;
}

size_t VtkType::size() const {
    switch (m_value) {
        case VtkType::Int64:
        case VtkType::UInt64:
        case VtkType::Float64:
            return 8;
        case VtkType::Int32:
        case VtkType::UInt32:
        case VtkType::Float32:
            return 4;
        case VtkType::Int16:
        case VtkType::UInt16:
            return 2;
        case VtkType::Int8:
        case VtkType::UInt8:
            return 1;
        default:
            return 0;
    }
}

VtkType VtkType::get_signed(size_t size) {
    switch (size) {
        case 1:  return VtkType::Int8;
        case 2:  return VtkType::Int16;
        case 4:  return VtkType::Int32;
        case 8:  return VtkType::Int64;
        default: return VtkType::Undefined;
    }
}

VtkType VtkType::get_unsigned(size_t size) {
    switch (size) {
        case 1:  return VtkType::UInt8;
        case 2:  return VtkType::UInt16;
        case 4:  return VtkType::UInt32;
        case 8:  return VtkType::UInt64;
        default: return VtkType::Undefined;
    }
}

VtkType VtkType::get_floating(size_t size) {
    switch (size) {
        case 4:  return VtkType::Float32;
        case 8:  return VtkType::Float64;
        default: return VtkType::Undefined;
    }
}

std::ostream& operator<<(std::ostream& os, const VtkType& vtk_type) {
    os << static_cast<const char *>(vtk_type);
    return os;
}

} // namespace zephyr::io