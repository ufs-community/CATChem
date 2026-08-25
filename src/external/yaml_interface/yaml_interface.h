#ifndef YAML_INTERFACE_H
#define YAML_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

// Node management functions
void* yaml_load_file(const char* filename);
void* yaml_load_string(const char* yaml_string);
void yaml_destroy_node(void* node_ptr);
void* yaml_sequence_to_map(void* node_ptr);

// Getter functions
bool yaml_get_string(void* node_ptr, const char* key, char* value, int max_len);
bool yaml_get_integer(void* node_ptr, const char* key, int* value);
bool yaml_get_real(void* node_ptr, const char* key, double* value);
bool yaml_get_logical(void* node_ptr, const char* key, bool* value);

// Array getter functions
bool yaml_get_real_array(void* node_ptr, const char* key, double* values, int max_size, int* actual_size);
bool yaml_get_integer_array(void* node_ptr, const char* key, int* values, int max_size, int* actual_size);
bool yaml_get_string_array(void* node_ptr, const char* key, char* values, int max_strings, int max_len,
                           int* actual_size);

// Utility functions
bool yaml_has_key(void* node_ptr, const char* key);
int yaml_get_size(void* node_ptr);
bool yaml_is_sequence(void* node_ptr);
bool yaml_is_map(void* node_ptr);
bool yaml_get_all_keys(void* node_ptr, char* keys, int max_keys, int max_key_len, int* actual_count);

// Setter functions
bool yaml_set_string(void* node_ptr, const char* key, const char* value);
bool yaml_set_integer(void* node_ptr, const char* key, int value);
bool yaml_set_real(void* node_ptr, const char* key, double value);
bool yaml_set_logical(void* node_ptr, const char* key, bool value);

// File operations
bool yaml_save_file(void* node_ptr, const char* filename);

#ifdef __cplusplus
}
#endif

#endif // YAML_INTERFACE_H
