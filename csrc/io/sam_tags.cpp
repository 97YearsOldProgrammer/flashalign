#include "sam_tags.h"

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace fa { namespace cpu { namespace io {

namespace {

bool tag_name_ok(std::string_view field) {
    const bool a = (field[0] >= 'A' && field[0] <= 'Z') ||
                   (field[0] >= 'a' && field[0] <= 'z');
    const bool b = (field[1] >= 'A' && field[1] <= 'Z') ||
                   (field[1] >= 'a' && field[1] <= 'z') ||
                   (field[1] >= '0' && field[1] <= '9');
    return a && b;
}

// Strict decimal integer (optional sign, digits only) within [lo, hi].
bool parse_int(std::string_view v, long long lo, long long hi, long long& out) {
    if (v.empty()) return false;
    size_t i = 0;
    bool neg = false;
    if (v[0] == '-' || v[0] == '+') { neg = (v[0] == '-'); i = 1; }
    if (i == v.size()) return false;
    long long value = 0;
    for (; i < v.size(); ++i) {
        if (v[i] < '0' || v[i] > '9') return false;
        const int digit = v[i] - '0';
        // Reject overflow rather than wrap.
        if (value > (4294967295LL - digit) / 10) return false;
        value = value * 10 + digit;
    }
    if (neg) value = -value;
    if (value < lo || value > hi) return false;
    out = value;
    return true;
}

bool parse_float(std::string_view v, float& out) {
    if (v.empty()) return false;
    // strtof skips leading whitespace and accepts it as a number; SAM does not.
    if (!((v[0] >= '0' && v[0] <= '9') || v[0] == '-' || v[0] == '+' ||
          v[0] == '.')) {
        return false;
    }
    const std::string copy(v);
    char* end = nullptr;
    const float value = std::strtof(copy.c_str(), &end);
    if (end != copy.c_str() + copy.size()) return false;
    out = value;
    return true;
}

void put_le(std::string& out, unsigned long long value, int width) {
    for (int k = 0; k < width; ++k)
        out.push_back(static_cast<char>((value >> (8 * k)) & 0xff));
}

void put_float(std::string& out, float value) {
    char raw[sizeof(float)];
    std::memcpy(raw, &value, sizeof(float));
    out.append(raw, sizeof(float));
}

// The smallest BAM integer code the value fits, as htslib chooses.
void put_int(std::string& out, long long value) {
    char code;
    int width;
    if (value >= 0) {
        if (value <= 255)        { code = 'C'; width = 1; }
        else if (value <= 65535) { code = 'S'; width = 2; }
        else                     { code = 'I'; width = 4; }
    } else {
        if (value >= -128)        { code = 'c'; width = 1; }
        else if (value >= -32768) { code = 's'; width = 2; }
        else                      { code = 'i'; width = 4; }
    }
    out.push_back(code);
    put_le(out, static_cast<unsigned long long>(value), width);
}

bool array_element_bounds(char subtype, long long& lo, long long& hi) {
    switch (subtype) {
        case 'c': lo = -128; hi = 127; return true;
        case 'C': lo = 0; hi = 255; return true;
        case 's': lo = -32768; hi = 32767; return true;
        case 'S': lo = 0; hi = 65535; return true;
        case 'i': lo = -2147483648LL; hi = 2147483647LL; return true;
        case 'I': lo = 0; hi = 4294967295LL; return true;
        default: return false;
    }
}

// Validates one tag field and, when `out` is non-null, encodes it. One
// function for both, so validation and encoding cannot disagree.
bool convert_field(std::string_view field, std::string* out) {
    if (field.size() < 5 || field[2] != ':' || field[4] != ':') return false;
    if (!tag_name_ok(field)) return false;
    const char type = field[3];
    const std::string_view value = field.substr(5);
    std::string encoded;
    encoded.push_back(field[0]);
    encoded.push_back(field[1]);
    switch (type) {
        case 'A': {
            if (value.size() != 1) return false;
            if (value[0] < '!' || value[0] > '~') return false;
            encoded.push_back('A');
            encoded.push_back(value[0]);
            break;
        }
        case 'i': {
            long long parsed = 0;
            if (!parse_int(value, -2147483648LL, 4294967295LL, parsed))
                return false;
            put_int(encoded, parsed);
            break;
        }
        case 'f': {
            float parsed = 0.0f;
            if (!parse_float(value, parsed)) return false;
            encoded.push_back('f');
            put_float(encoded, parsed);
            break;
        }
        case 'Z': {
            for (char ch : value)
                if (ch < ' ' || ch > '~') return false;
            encoded.push_back('Z');
            encoded.append(value.data(), value.size());
            encoded.push_back('\0');
            break;
        }
        case 'H': {
            if (value.size() % 2 != 0) return false;
            for (char ch : value) {
                const bool hex = (ch >= '0' && ch <= '9') ||
                                 (ch >= 'A' && ch <= 'F') ||
                                 (ch >= 'a' && ch <= 'f');
                if (!hex) return false;
            }
            encoded.push_back('H');
            encoded.append(value.data(), value.size());
            encoded.push_back('\0');
            break;
        }
        case 'B': {
            if (value.empty()) return false;
            const char subtype = value[0];
            long long lo = 0, hi = 0;
            const bool integral = array_element_bounds(subtype, lo, hi);
            if (!integral && subtype != 'f') return false;
            const int width = (subtype == 'c' || subtype == 'C') ? 1
                            : (subtype == 's' || subtype == 'S') ? 2 : 4;
            // Parse every element before writing the count prefix.
            std::vector<long long> integers;
            std::vector<float> floats;
            std::string_view rest = value.substr(1);
            while (!rest.empty()) {
                if (rest[0] != ',') return false;
                rest.remove_prefix(1);
                const size_t comma = rest.find(',');
                const std::string_view element =
                    comma == std::string_view::npos ? rest
                                                    : rest.substr(0, comma);
                if (integral) {
                    long long parsed = 0;
                    if (!parse_int(element, lo, hi, parsed)) return false;
                    integers.push_back(parsed);
                } else {
                    float parsed = 0.0f;
                    if (!parse_float(element, parsed)) return false;
                    floats.push_back(parsed);
                }
                rest = comma == std::string_view::npos
                           ? std::string_view()
                           : rest.substr(comma);
            }
            encoded.push_back('B');
            encoded.push_back(subtype);
            const size_t count = integral ? integers.size() : floats.size();
            put_le(encoded, static_cast<unsigned long long>(count), 4);
            for (long long element : integers) {
                put_le(encoded, static_cast<unsigned long long>(element),
                       width);
            }
            for (float element : floats) put_float(encoded, element);
            break;
        }
        default: return false;
    }
    if (out) *out += encoded;
    return true;
}

bool convert_text(std::string_view text, std::string* out,
                  std::string* offending_field) {
    if (text.empty()) return true;
    size_t begin = 0;
    for (;;) {
        size_t end = text.find('\t', begin);
        if (end == std::string_view::npos) end = text.size();
        const std::string_view field = text.substr(begin, end - begin);
        if (!convert_field(field, out)) {
            if (offending_field) offending_field->assign(field);
            return false;
        }
        if (end == text.size()) break;
        begin = end + 1;
    }
    return true;
}

}  // namespace

bool sam_tag_text_valid(std::string_view text, std::string* offending_field) {
    return convert_text(text, nullptr, offending_field);
}

std::string sam_tag_text_to_bam(std::string_view text) {
    std::string out;
    std::string offending;
    if (!convert_text(text, &out, &offending)) {
        throw std::runtime_error("not SAM tag text: " + offending);
    }
    return out;
}

}}}  // namespace fa::cpu::io
