# Process Template Updates: Unified Configuration Architecture

## Overview

The process templates have been updated to implement a **unified configuration architecture** that properly integrates ProcessInterface with ConfigManager's hierarchical YAML system.

## Key Changes

### 1. **Template Files Updated**

- `process_interface.F90.j2` - Main process interface template
- `process_common.F90.j2` - Common types and configuration template

### 2. **New Architecture Components**

#### **Unified Process Configuration Type**
```fortran
type :: {{ config.class_name }}ProcessConfig
   character(len=64) :: process_name
   character(len=16) :: process_version
   logical :: is_active

   ! Process-specific configuration
   type({{ config.class_name }}Config) :: {{ config.name }}_config

   ! Scheme configurations
   type({{ config.class_name }}Scheme{{ scheme.class_name }}Config) :: {{ scheme.name }}_config

contains
   procedure :: load_from_config
   procedure :: validate
   procedure :: finalize
   procedure :: get_active_scheme_config
end type
```

#### **Proper parse_process_config Implementation**
```fortran
subroutine parse_{{ config.name }}_config(this, config_data, error_manager)
   ! Uses unified configuration loader from Common module
   call this%process_config%load_from_config(config_data, error_manager)
end subroutine
```

#### **Hierarchical YAML Support**
```fortran
subroutine {{ config.name }}_process_load_config(this, config_data, error_handler)
   ! Parses processes.{{ config.name }}.* structure
   ! Loads process metadata, main config, and scheme-specific configs
end subroutine
```

### 3. **Template Generation Flow**

1. **ConfigManager** loads hierarchical YAML
2. **ProcessInterface.parse_process_config()** delegates to unified loader
3. **{{ config.class_name }}ProcessConfig.load_from_config()** parses the YAML structure
4. **Scheme methods** access config via `this%process_config%{{ scheme.name }}_config`

## YAML Configuration Structure

The templates now expect this hierarchical structure:

```yaml
processes:
  {{ config.name }}:
    # Process metadata
    name: "{{ config.name }}"
    version: "1.0.0"
    active: true

    # Main process configuration
    scheme: "{{ scheme.name }}"
    dt_min: 60.0
    dt_max: 3600.0

    # Species configuration
    species:
      - name: "species1"
        active: true

    # Scheme-specific configuration (optional)
    scheme_config:
      scale_factor: 1.0
      # scheme-specific parameters

    # Diagnostic configuration
    diagnostics:
      diagnostic_name:
        active: true
        units: "units"
```

## Benefits

### ✅ **Problems Solved**

1. **Missing parse_process_config** - Now properly implemented in templates
2. **ConfigManager Integration** - Templates work with hierarchical YAML structure
3. **Template Consistency** - All generated processes follow unified architecture
4. **Configuration Validation** - Proper error handling and validation

### 🎯 **Architecture Advantages**

1. **Clean Separation** - ConfigManager ↔ Unified Config ↔ Process-specific types
2. **Flexible Configuration** - Support for process-level and scheme-specific settings
3. **Extensible Design** - Easy to add new processes and schemes
4. **Consistent Interface** - All processes follow the same configuration pattern

## Usage

### **Regenerating Processes**

When using the process generator with these updated templates:

1. **Existing YAML configs** need minor updates to match expected structure
2. **Generated code** will include proper ConfigManager integration
3. **ProcessInterface** will have working `parse_process_config` implementation

### **Configuration Loading**

The generated code will:

1. Parse `processes.{{ config.name }}.*` from main YAML
2. Load process metadata (name, version, active status)
3. Initialize process-specific configuration
4. Load scheme-specific parameters from `scheme_config` section
5. Validate all configuration with StateManager

## Example

For the seasalt process:

```fortran
! Generated interface includes:
type(SeaSaltProcessConfig) :: process_config

! Configuration loading:
call this%process_config%load_from_config(config_data, error_manager)

! Scheme access:
select case (trim(this%process_config%seasalt_config%scheme))
case ('gong97')
   call compute_gong97(..., this%process_config%gong97_config, ...)
end select
```

## Migration

### **For Existing Processes**

1. **Update YAML structure** to match new hierarchy
2. **Regenerate process files** with updated templates
3. **Update scheme calls** to use unified configuration access
4. **Test configuration loading** with ConfigManager

### **For New Processes**

1. **Use updated templates** from process generator
2. **Follow YAML structure** shown in examples
3. **Implement scheme modules** with expected interfaces
4. **All configuration handling** is automatically generated

This unified architecture provides a robust, extensible foundation for CATChem process configuration management!
