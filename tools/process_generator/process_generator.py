#!/usr/bin/env python3
"""
CATChem Process Generator Tool

A comprehensive tool for generating standardized process implementations in CATChem.
Uses YAML configuration files and Jinja2 templates to create complete process
modules, schemes, documentation, tests, and CMake integration.

This tool follows the Process Infrastructure Guide and creates processes
that are compatible with the modern CATChem architecture.

Usage:
    python process_generator.py generate --config my_process.yaml
    python process_generator.py validate --config my_process.yaml
    python process_generator.py template --type process --output template.yaml

Author: CATChem Development Team
License: Apache 2.0
"""

import argparse
import dataclasses
import logging
import os
import sys
import yaml
from pathlib import Path
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
from jinja2 import Environment, FileSystemLoader, select_autoescape
import json
from datetime import datetime
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('ProcessGenerator')


def _extract_field_names(field_list) -> List[str]:
    """Extract string field names from a list of strings or dicts and normalize to standard MetState names."""
    alias_map = {
        'LEAF_AREA_INDEX': 'LAI',
        'SOLAR_ZENITH_ANGLE': 'SUNCOS',
        'TEMPERATURE': 'T',
        'PRESSURE': 'PMID',
        'SURFACE_PRESSURE': 'PS',
        'HUMIDITY': 'QV',
        'AIR_DENSITY': 'AIRDEN',
        'WIND_SPEED': 'U10M',
        'FRICTION_VELOCITY': 'USTAR',
    }
    names = []
    seen = set()
    if not field_list:
        return names
    for item in field_list:
        name = None
        if isinstance(item, dict):
            name = item.get('name') or item.get('field') or item.get('variable_name')
        elif isinstance(item, str):
            name = item
        else:
            name = str(item)
        if name:
            norm_name = str(name).strip().upper()
            norm_name = alias_map.get(norm_name, norm_name)
            if norm_name not in seen:
                seen.add(norm_name)
                names.append(norm_name)
    return names


@dataclass
class ProcessBehavior:
    """Configuration for process behavior patterns."""
    type: str = "source"                    # source, sink, transformation, transport
    tendency_mode: str = "additive"         # additive, replacement, multiplicative
    species_filter: Dict[str, Any] = field(default_factory=dict)
    tendency_calculation: str = "rates"     # rates, concentrations, deltas
    timestep_dependency: str = "independent"  # independent, dependent, adaptive
    spatial_scope: str = "column"           # column, global, regional
    parallelization: str = "column"         # column, species, domain
    memory_requirements: str = "low"        # low, medium, high
    gas_aero_differentiation: bool = False  # Whether to separate gas/aero schemes


@dataclass
class SchemeBehavior:
    """Configuration for scheme-specific behavior."""
    output_format: str = "rates_2d"         # rates_2d, concentrations_1d, delta_concentrations
    input_requirements: List[str] = field(default_factory=list)  # vertical_profile, surface_properties, etc.


@dataclass
class MetFieldClassification:
    """Classification for meteorological field types based on MetState definition."""

    def __init__(self, metstate_file: Optional[str] = None):
        """Initialize with optional MetState file path for automatic field discovery."""
        self.metstate_file = metstate_file
        self._fields_cache = None

        # Fallback hardcoded lists for when MetState file is not available
        self._fallback_2d_surface = [
            'FROCEAN', 'FRSEAICE', 'FRLAKE', 'FRLAND', 'SST', 'TSK', 'SKINTEMP',
            'U10M', 'V10M', 'T2M', 'Q2M', 'PS', 'SLP', 'USTAR', 'PBLH',
            'SNOWH', 'ALBEDO', 'EMISS', 'HFX', 'QFX', 'LH', 'MOL', 'TS', 'TSKIN',
            'QV2M', 'PHIS', 'SUNCOS', 'SWGDN', 'HFLUX', 'EFLUX', 'PRECCON',
            'PRECLSC', 'PRECANV', 'TO3', 'TROPP', 'TropHt', 'Z0', 'LAI', 'AREA_M2'
        ]
        self._fallback_3d_atmospheric = [
            'T', 'QV', 'P', 'PLE', 'DELP', 'U', 'V', 'OMEGA', 'W',
            'RH', 'CLOUD', 'QC', 'QR', 'QI', 'QS', 'QG', 'PMID', 'AIRDEN',
            'THETA', 'SPHU', 'AIRNUMDEN', 'MAIRDEN', 'AVGW', 'DELP_DRY',
            'DAIRMASS', 'AIRVOL', 'PEDGE_DRY', 'PEDGE', 'PMID_DRY', 'Z',
            'ZMID', 'BXHEIGHT', 'TV', 'CLDF', 'CMFMC', 'DQRCU', 'DQRLSAN',
            'DTRAIN', 'QL', 'PFICU', 'PFILSAN', 'PFLCU', 'PFLLSAN',
            'TAUCLI', 'TAUCLW', 'F_OF_PBL', 'F_UNDER_PBLTOP'
        ]
        self._fallback_categorical = [
            'SOILM', 'FRLANDUSE', 'FRSOIL', 'FRLAI', 'FRZ0', 'ILAND', 'LU_INDEX',
            'VEGFRA', 'SOILTYP', 'SLOPETYP', 'XLAND', 'IVGTYP', 'ISLTYP',
            'VEGETATION_TYPE'
        ]

    def _parse_metstate_fields(self) -> Dict[str, tuple]:
        """Parse MetState file to extract field definitions."""
        if self._fields_cache is not None:
            return self._fields_cache

        if not self.metstate_file or not Path(self.metstate_file).exists():
            logger.warning(f"MetState file not found: {self.metstate_file}, using fallback classifications")
            self._fields_cache = {}
            return self._fields_cache

        try:
            # Use the same parsing logic as generate_metstate_macros.py
            fields_cache = {}

            with open(self.metstate_file, 'r') as f:
                lines = f.readlines()

            in_type = False
            for line in lines:
                if 'TYPE, PUBLIC :: MetStateType' in line:
                    in_type = True
                elif in_type and 'end type' in line.lower():
                    break
                elif in_type:
                    # Match real(fp), allocatable :: name(dimensions)
                    import re
                    m = re.match(r'\s*REAL\(fp\),\s*ALLOCATABLE\s*::\s*(\w+)\s*\(([^)]*)\)', line, re.IGNORECASE)
                    if m:
                        name = m.group(1)
                        dims = m.group(2)
                        rank = dims.count(',') + 1
                        # Check if this is an edge field (nz+1 dimension)
                        is_edge = 'nz+1' in line or 'nlevs+1' in line.lower()
                        fields_cache[name] = (rank, dims, is_edge)

            self._fields_cache = fields_cache
            logger.info(f"Parsed {len(fields_cache)} fields from MetState file")
            return self._fields_cache

        except Exception as e:
            logger.warning(f"Error parsing MetState file {self.metstate_file}: {e}, using fallback classifications")
            self._fields_cache = {}
            return self._fields_cache

    def get_field_type(self, field_name: str) -> str:
        """Get the type of a meteorological field."""

        # Handle special non-meteorological fields first
        if field_name == 'TSTEP':
            return 'special_timestep'

        fields = self._parse_metstate_fields()

        if field_name in fields:
            rank, dims, is_edge = fields[field_name]

            # Categorize based on rank and known categorical fields
            if rank == 2:
                return '2d_surface'
            elif rank == 3:
                # Check if it's an edge field
                if is_edge:
                    return '3d_edge'
                # Check if it's a categorical 3D field (special dimensions)
                elif field_name in ['SOILM', 'FRLANDUSE', 'FRSOIL', 'FRLAI', 'FRZ0', 'ILAND']:
                    return 'categorical'
                else:
                    return '3d_atmospheric'
            else:
                # Scalar or other rank
                return '2d_surface'  # Default
        else:
            # Use fallback classification
            if field_name in self._fallback_2d_surface:
                return '2d_surface'
            elif field_name in self._fallback_3d_atmospheric:
                return '3d_atmospheric'
            elif field_name in self._fallback_categorical:
                return 'categorical'
            else:
                # Default to 2d_surface for unknown fields
                return '2d_surface'

    def get_array_size(self, field_name: str, affects_full_column: bool) -> str:
        """Get the array size specification for a field based on affects_full_column."""
        field_type = self.get_field_type(field_name)

        if field_type == 'categorical':
            return '(:)'  # Always 1D array for categorical
        elif field_type == '2d_surface':
            return ''  # Always scalar for 2D surface fields
        elif field_type == '3d_edge':
            if affects_full_column:
                return '(:)'  # 1D array for full column (nz+1 size)
            else:
                return ''  # Scalar for surface only
        elif field_type == '3d_atmospheric':
            if affects_full_column:
                return '(:)'  # 1D array for full column
            else:
                return ''  # Scalar for surface only
        else:
            return ''  # Default to scalar

    def get_all_2d_fields(self) -> List[str]:
        """Get all 2D surface fields."""
        fields = self._parse_metstate_fields()
        result = []

        # From parsed fields
        for field_name, field_info in fields.items():
            rank = field_info[0]  # Handle both 2-tuple and 3-tuple formats
            if rank == 2:
                result.append(field_name)

        # Add fallback fields not found in parsed file
        for field_name in self._fallback_2d_surface:
            if field_name not in result:
                result.append(field_name)

        return sorted(result)

    def get_all_3d_atmospheric_fields(self) -> List[str]:
        """Get all 3D atmospheric fields (excluding categorical and edge)."""
        fields = self._parse_metstate_fields()
        result = []

        # From parsed fields
        categorical_3d = {'SOILM', 'FRLANDUSE', 'FRSOIL', 'FRLAI', 'FRZ0'}
        for field_name, field_info in fields.items():
            rank = field_info[0]  # Handle both 2-tuple and 3-tuple formats
            is_edge = field_info[2] if len(field_info) >= 3 else False
            if rank == 3 and field_name not in categorical_3d and not is_edge:
                result.append(field_name)

        # Add fallback fields not found in parsed file
        for field_name in self._fallback_3d_atmospheric:
            if field_name not in result:
                result.append(field_name)

        return sorted(result)

    def get_all_3d_edge_fields(self) -> List[str]:
        """Get all 3D edge fields (nz+1 dimension)."""
        fields = self._parse_metstate_fields()
        result = []

        # From parsed fields - only those with is_edge=True
        for field_name, field_info in fields.items():
            rank = field_info[0]
            is_edge = field_info[2] if len(field_info) >= 3 else False
            if rank == 3 and is_edge:
                result.append(field_name)

        return sorted(result)

    def get_data_type(self, field_name: str) -> str:
        """Get the Fortran data type for a meteorological field."""
        # Boolean/logical fields - fields that represent true/false values
        boolean_fields = {
            'IsSnow', 'IsIce', 'IsLand', 'ISICE', 'ISSNOW', 'ISLAND',
            'is_snow', 'is_ice', 'is_land'
        }

        # Character/string fields - fields that represent text/names
        character_fields = {
            'LUCNAME', 'State', 'lucname', 'state'
        }

        # Integer fields - fields that represent integer values
        integer_fields = {
            'ILAND', 'iland', 'LWI', 'lwi', 'DLUSE', 'dluse',
            'DSOILTYPE', 'dsoiltype', 'nLNDTYPE', 'nlndtype',
            'TropLev', 'troplev'
        }

        if field_name in boolean_fields:
            return 'logical'
        elif field_name in character_fields:
            return 'character(len=255)'
        elif field_name in integer_fields:
            return 'integer'
        else:
            return 'real(fp)'

    def get_species_property_data_type(self, property_name: str) -> str:
        """Get the Fortran data type for a species property."""
        # Boolean/logical species properties - properties that represent true/false values
        boolean_properties = {
            'is_dust', 'is_seasalt', 'is_gas', 'is_aerosol', 'is_tracer',
            'is_transported', 'is_wet_scavenged', 'is_dry_deposited',
            'wd_LiqAndGas'  # Add wet deposition logical property
        }

        # Character/string species properties - properties that represent text/names
        string_properties = {
            'short_name', 'name', 'long_name', 'formula', 'chem_formula'
        }

        if property_name in boolean_properties:
            return 'logical'
        elif property_name in string_properties:
            return 'character(len=32)'
        else:
            return 'real(fp)'

    def get_species_property_dimensions(self, property_name: str) -> str:
        """Get the Fortran array dimensions for a species property."""
        # Properties with special dimensions
        if property_name == 'wd_rainouteff':
            return '(:,:)'  # 2D array: (n_species, 3)
        else:
            return '(:)'    # Default: 1D array (n_species)

    def get_species_property_allocation_size(self, property_name: str) -> str:
        """Get the Fortran allocation size for a species property."""
        # Properties with special dimensions
        if property_name == 'wd_rainouteff':
            return '(this%{{ config.name }}_config%n_species, 3)'
        else:
            return '(this%{{ config.name }}_config%n_species)'

    def get_all_categorical_fields(self) -> List[str]:
        """Get all categorical fields."""
        fields = self._parse_metstate_fields()
        result = []

        # From parsed fields
        categorical_3d = {'SOILM', 'FRLANDUSE', 'FRSOIL', 'FRLAI', 'FRZ0'}
        for field_name, field_info in fields.items():
            rank = field_info[0]  # Handle both 2-tuple and 3-tuple formats
            if field_name in categorical_3d:
                result.append(field_name)

        # Add fallback fields not found in parsed file
        for field_name in self._fallback_categorical:
            if field_name not in result:
                result.append(field_name)

        return sorted(result)


@dataclass
class SchemeConfig:
    """Configuration for a process scheme."""
    name: str
    class_name: str
    description: str
    author: str = ""
    reference: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)
    required_met_fields: List[str] = field(default_factory=list)
    required_species_properties: List[str] = field(default_factory=list)
    required_constants: List[str] = field(default_factory=list)
    required_time_parameters: List[str] = field(default_factory=list)
    scheme_diagnostics: List[Dict[str, str]] = field(default_factory=list)
    persistent_state_variables: List[Dict[str, Any]] = field(default_factory=list)
    algorithm_type: str = "explicit"
    affects_full_column: bool = False  # Whether scheme affects full atmospheric column
    scheme_type: str = ""  # Optional legacy field
    scheme_behavior: Optional[SchemeBehavior] = None
    gas_or_aero: str = "both"  # New field: gas, aero, or both


@dataclass
class ProcessConfig:
    """Main process configuration."""
    name: str
    description: str
    class_name: str
    author: str
    version: str = "1.0.0"
    license: str = "Apache 2.0"

    # Process behavior configuration (replaces hardcoded process_type)
    process_behavior: Optional[ProcessBehavior] = None

    # Legacy field for backward compatibility
    process_type: str = "generic"

    is_multiphase: bool = False
    has_size_bins: bool = False
    supports_vectorization: bool = True
    species: List[str] = field(default_factory=list)
    size_bins: Optional[Dict[str, Any]] = None
    phases: List[str] = field(default_factory=lambda: ['gas'])
    schemes: List[SchemeConfig] = field(default_factory=list)
    default_scheme: str = ""
    required_met_fields: List[str] = field(default_factory=list)
    optional_met_fields: List[str] = field(default_factory=list)
    required_constants: List[str] = field(default_factory=list)
    required_chem_fields: List[str] = field(default_factory=list)
    diagnostics: List[Dict[str, str]] = field(default_factory=list)
    diagnostic_species: List[str] = field(default_factory=lambda: ["All"])  # Default to all species for diagnostics
    timestep_dependency: str = "independent"
    parallelization: str = "column"
    memory_requirements: str = "low"
    generate_tests: bool = True
    generate_docs: bool = True
    generate_examples: bool = False
    output_dir: str = ""
    src_base_dir: str = "src/process"

    @property
    def gas_aero_differentiation(self) -> bool:
        """Determine if gas/aero differentiation is enabled from process_behavior."""
        if self.process_behavior and hasattr(self.process_behavior, 'gas_aero_differentiation'):
            return self.process_behavior.gas_aero_differentiation
        return False

    @property
    def enable_column_processing(self) -> bool:
        """Determine if column processing is enabled based on parallelization strategy."""
        if self.process_behavior and self.process_behavior.parallelization:
            return self.process_behavior.parallelization == "column"
        return True  # Default to column processing if not specified


class ProcessValidationError(Exception):
    """Exception raised for process configuration validation errors."""
    pass


class ProcessGenerator:
    """Main process generator class."""

    def __init__(self, template_dir: Optional[str] = None, metstate_file: Optional[str] = None):
        """Initialize the process generator.

        Args:
            template_dir: Directory containing Jinja2 templates. If None,
                         uses default templates in same directory as this script.
            metstate_file: Path to MetState file for automatic field discovery. If None,
                          tries to find it automatically relative to the script location.
        """
        if template_dir is None:
            template_dir = str(Path(__file__).parent / "templates")

        self.template_dir = Path(template_dir)

        # Try to find MetState file automatically if not provided
        if metstate_file is None:
            script_dir = Path(__file__).resolve().parent
            # Look for metstate_mod.F90 relative to process generator
            potential_paths = [
                script_dir.parent.parent / "src" / "core" / "metstate_mod.F90",  # From tools/process_generator
                script_dir / "../../src/core/metstate_mod.F90",  # Alternative relative path
            ]
            for path in potential_paths:
                if path.exists():
                    metstate_file = str(path)
                    logger.info(f"Found MetState file automatically: {metstate_file}")
                    break

        self.metstate_file = metstate_file
        if metstate_file and not Path(metstate_file).exists():
            logger.warning(f"MetState file not found: {metstate_file}")

        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=select_autoescape(['html', 'xml']),
            trim_blocks=True,
            lstrip_blocks=True
        )

        # Add custom filters
        self.env.filters['upper_snake'] = self._upper_snake_case
        self.env.filters['lower_snake'] = self._lower_snake_case
        self.env.filters['pascal_case'] = self._pascal_case
        self.env.filters['camel_case'] = self._camel_case
        self.env.filters['fortran_string'] = self._fortran_string
        self.env.filters['fortran_boolean'] = self._fortran_boolean
        self.env.filters['infer_diagnostic_type'] = self._infer_diagnostic_type
        self.env.filters['infer_diagnostic_properties'] = self._infer_diagnostic_properties
        self.env.filters['analyze_required_dimensions'] = self._analyze_required_dimensions
        self.env.filters['fortran_array_constructor'] = self._fortran_array_constructor

        # Add custom tests
        self.env.tests['list_type'] = lambda val: isinstance(val, list)

        # Add a filter to get all required met fields for a scheme
        def get_all_met_fields_filter(scheme):
            """Jinja2 filter to get all required meteorological fields for a scheme."""
            context = self.env.globals.get('config', {})
            all_fields = set()

            if hasattr(context, 'required_met_fields') and context.required_met_fields:
                all_fields.update(_extract_field_names(context.required_met_fields))

            scheme_fields = None
            if isinstance(scheme, dict):
                scheme_fields = scheme.get('required_met_fields')
            elif hasattr(scheme, 'required_met_fields'):
                scheme_fields = scheme.required_met_fields

            if scheme_fields:
                all_fields.update(_extract_field_names(scheme_fields))

            return sorted(list(all_fields))

        # Add filter for scheme-only met fields (excludes process-level fields)
        def get_scheme_only_met_fields_filter(scheme, context=None):
            """Filter for scheme-only meteorological fields (excludes process-level fields)."""
            scheme_fields = set()

            if isinstance(scheme, dict):
                if 'required_met_fields' in scheme and scheme['required_met_fields']:
                    scheme_fields.update(_extract_field_names(scheme['required_met_fields']))
            elif hasattr(scheme, 'required_met_fields') and scheme.required_met_fields:
                scheme_fields.update(_extract_field_names(scheme.required_met_fields))

            return sorted(list(scheme_fields))

        self.env.filters['all_required_met_fields'] = get_all_met_fields_filter
        self.env.filters['scheme_only_met_fields'] = get_scheme_only_met_fields_filter

    @staticmethod
    def _upper_snake_case(s: str) -> str:
        """Convert string to UPPER_SNAKE_CASE."""
        return s.upper().replace(' ', '_').replace('-', '_')

    @staticmethod
    def _lower_snake_case(s: str) -> str:
        """Convert string to lower_snake_case."""
        return s.lower().replace(' ', '_').replace('-', '_')

    @staticmethod
    def _pascal_case(s: str) -> str:
        """Convert string to PascalCase."""
        return ''.join(word.capitalize() for word in s.replace('_', ' ').replace('-', ' ').split())

    @staticmethod
    def _camel_case(s: str) -> str:
        """Convert string to camelCase."""
        pascal = ProcessGenerator._pascal_case(s)
        return pascal[0].lower() + pascal[1:] if pascal else ""

    @staticmethod
    def _fortran_string(s: str, length: int = 64) -> str:
        """Format string for Fortran character declaration with proper padding."""
        # Pad or truncate string to exact length for array constructors
        padded = s[:length].ljust(length)
        return f"'{padded}'"

    @staticmethod
    def _fortran_boolean(b: bool) -> str:
        """Convert boolean to Fortran logical."""
        return ".true." if b else ".false."

    @staticmethod
    def _fortran_array_constructor(values, suffix: str = '_fp') -> str:
        """Convert a Python list to a Fortran array constructor string.

        Examples:
            [1.0, 2.0, 3.0] -> '(/ 1.0_fp, 2.0_fp, 3.0_fp /)'
            [1, 2, 3]       -> '(/ 1, 2, 3 /)'
        """
        if not isinstance(values, list) or len(values) == 0:
            return '0.0' + suffix
        # Detect element type from first element
        if all(isinstance(v, bool) for v in values):
            items = ['.true.' if v else '.false.' for v in values]
        elif all(isinstance(v, int) and not isinstance(v, bool) for v in values):
            items = [str(v) for v in values]
        else:
            # Treat as real
            items = [str(v) + suffix for v in values]
        return '(/ ' + ', '.join(items) + ' /)'

    @staticmethod
    def _infer_state_variable_type(default_value: Any, name: str) -> str:
        """Infer Fortran type from default value and variable name."""
        if isinstance(default_value, bool):
            return "logical"
        elif isinstance(default_value, int):
            return "integer"
        elif isinstance(default_value, (float, int)):
            return "real(fp)"
        elif name.endswith('(:)') or '(:)' in name:
            # Array variable - infer base type from default
            if isinstance(default_value, bool):
                return "logical"
            elif isinstance(default_value, int):
                return "integer"
            else:
                return "real(fp)"
        else:
            return "real(fp)"  # Default to real

    @staticmethod
    def _get_state_variable_dimensions(name: str) -> str:
        """Get array dimensions from variable name."""
        if '(:)' in name:
            # Extract dimensions - for now support (:) which means allocatable 1D
            return "(:)"
        else:
            return ""  # Scalar

    @staticmethod
    def _clean_state_variable_name(name: str) -> str:
        """Clean variable name by removing dimension specifications."""
        return name.replace('(:)', '').strip()

    def _infer_diagnostic_type(self, diagnostic: Dict[str, Any], config: ProcessConfig, scheme_config: SchemeConfig = None) -> str:
        """Infer diagnostic data type from configuration and context."""
        result = self._infer_diagnostic_properties(diagnostic, config, scheme_config)
        return result['data_type']

    def _infer_diagnostic_properties(self, diagnostic: Dict[str, Any], config: ProcessConfig, scheme_config: SchemeConfig = None) -> Dict[str, Any]:
        """Infer diagnostic data type and dimensions from configuration and context."""
        name = diagnostic.get('name', '')
        units = diagnostic.get('units', '')
        description = diagnostic.get('description', '')

        # Default result structure
        result = {
            'data_type': 'DIAG_REAL_2D',
            'dimensions': ['nx', 'ny'],
            'dimension_vars': ['dims_2d'],
            'fortran_dims': 'dims_2d',
            'dimension_source': 'grid_manager',  # How to get the dimensions
            'dimension_type': 'scalar',  # 'scalar', '1d', '2d', '3d'
            'dimension_name': None       # Primary dimension name for 1D arrays
        }

        # 1. Check for explicit specifications
        if 'data_type' in diagnostic:
            result['data_type'] = diagnostic['data_type']
        if 'dimensions' in diagnostic:
            result['dimensions'] = diagnostic['dimensions']
            result['fortran_dims'] = self._format_fortran_dims(diagnostic['dimensions'])
            # Update dimension type based on number of dimensions
            if len(diagnostic['dimensions']) == 0:
                result['dimension_type'] = 'scalar'
            elif len(diagnostic['dimensions']) == 1:
                result['dimension_type'] = '1d'
                result['dimension_name'] = diagnostic['dimensions'][0]
            elif len(diagnostic['dimensions']) == 2:
                result['dimension_type'] = '2d'
            else:
                result['dimension_type'] = '3d'
            return result

        # 2. Infer from field name and description patterns (continuous variables)
        name_lower = name.lower()
        desc_lower = description.lower()

        # Check for combined level AND species patterns for 4D diagnostics
        has_level_pattern = ('_per_level' in name or '_profile' in name or '_vertical' in name or
                           '_column' in name or '_layer' in name or
                           'level' in desc_lower or 'levels' in desc_lower or 'vertical' in desc_lower or
                           'profile' in desc_lower or 'column' in desc_lower or 'layer' in desc_lower or
                           'atmospheric' in desc_lower)

        has_species_pattern = ('_per_bin' in name or '_per_species' in name or '_per_mode' in name or '_distribution' in name or
                             'per bin' in desc_lower or 'per species' in desc_lower or 'per mode' in desc_lower or
                             'distribution' in desc_lower or 'size resolved' in desc_lower)

        # Priority 1: Check for combined level AND species patterns for 3D level diagnostics
        if has_level_pattern and has_species_pattern:
            result.update({
                'data_type': 'DIAG_REAL_3D',
                'dimensions': ['nx', 'ny', 'nz'],
                'dimension_vars': ['dims_3d_levels'],
                'fortran_dims': 'dims_3d_levels',
                'dimension_source': 'grid_manager',
                'dimension_type': '3d_levels_species',
                'dimension_name': 'levels_with_species'
            })

        # Priority 2: Check for species/bin/distribution patterns only
        elif has_species_pattern:
            result.update({
                'data_type': 'DIAG_REAL_3D',
                'dimensions': ['nx', 'ny', 'n_species'],
                'dimension_vars': ['dims_3d_species'],
                'fortran_dims': 'dims_3d_species',
                'dimension_source': 'process_config',
                'dimension_type': '1d',
                'dimension_name': 'n_species'
            })

        # Priority 3: Check for level/vertical patterns only
        elif has_level_pattern:
            result.update({
                'data_type': 'DIAG_REAL_3D',
                'dimensions': ['nx', 'ny', 'nz'],
                'dimension_vars': ['dims_3d_levels'],
                'fortran_dims': 'dims_3d_levels',
                'dimension_source': 'grid_manager',
                'dimension_type': '1d',
                'dimension_name': 'n_levels'
            })

        elif '_per_soil_layer' in name or '_soil_profile' in name:
            result.update({
                'data_type': 'DIAG_REAL_3D',
                'dimensions': ['nx', 'ny', 'n_soil_layers'],
                'dimension_vars': ['dims_soil'],
                'fortran_dims': 'dims_soil',
                'dimension_source': 'process_config',
                'dimension_type': '1d',
                'dimension_name': 'n_soil_layers'
            })

        # Check for column integrated patterns in name or description
        elif (('_column_integrated' in name or '_vertically_integrated' in name) or
              ('column integrated' in desc_lower or 'vertically integrated' in desc_lower or
               'integrated' in desc_lower)):
            # Column-integrated quantity - still 2D but from 3D process
            result.update({
                'data_type': 'DIAG_REAL_2D',
                'dimensions': ['nx', 'ny'],
                'dimension_vars': ['dims_2d'],
                'fortran_dims': 'dims_2d',
                'dimension_source': 'grid_manager',
                'dimension_type': 'scalar',
                'dimension_name': None
            })

        # Check for flux/total patterns in name or description
        elif (('_total' in name or '_integrated' in name or 'flux' in name) or
              ('total' in desc_lower or 'integrated' in desc_lower or 'flux' in desc_lower or
               'surface' in desc_lower or 'emission' in desc_lower)):
            # Check if it's a surface flux or column-integrated quantity
            if scheme_config and getattr(scheme_config, 'affects_full_column', False):
                # 3D process producing 2D output (column-integrated)
                result.update({
                    'data_type': 'DIAG_REAL_2D',
                    'dimensions': ['nx', 'ny'],
                    'dimension_vars': ['dims_2d'],
                    'fortran_dims': 'dims_2d',
                    'dimension_source': 'grid_manager',
                    'dimension_type': 'scalar',
                    'dimension_name': None
                })
            else:
                # Surface process producing 2D output
                result.update({
                    'data_type': 'DIAG_REAL_2D',
                    'dimensions': ['nx', 'ny'],
                    'dimension_vars': ['dims_2d'],
                    'fortran_dims': 'dims_2d',
                    'dimension_source': 'grid_manager',
                    'dimension_type': 'scalar',
                    'dimension_name': None
                })

        # 4. Infer from units
        elif units in ['unitless', 'dimensionless', 'index', 'category', '1', '-']:
            # Likely categorical or index data
            if 'probability' in name.lower() or 'fraction' in name.lower():
                result.update({
                    'data_type': 'DIAG_REAL_2D',
                    'dimensions': ['nx', 'ny'],
                    'dimension_vars': ['dims_2d'],
                    'fortran_dims': 'dims_2d',
                    'dimension_source': 'grid_manager',
                    'dimension_type': 'scalar',
                    'dimension_name': None
                })
            else:
                result.update({
                    'data_type': 'DIAG_INTEGER_2D',
                    'dimensions': ['nx', 'ny'],
                    'dimension_vars': ['dims_2d'],
                    'fortran_dims': 'dims_2d',
                    'dimension_source': 'grid_manager',
                    'dimension_type': 'scalar',
                    'dimension_name': None
                })

        # 5. Infer from process characteristics
        elif hasattr(config, 'process_type'):
            if config.process_type in ['emission', 'deposition']:
                # Most emissions/deposition are surface processes
                if scheme_config and getattr(scheme_config, 'affects_full_column', False):
                    # Full column emissions (rare)
                    result.update({
                        'data_type': 'DIAG_REAL_3D',
                        'dimensions': ['nx', 'ny', 'nz'],
                        'dimension_vars': ['dims_3d_levels'],
                        'fortran_dims': 'dims_3d_levels',
                        'dimension_source': 'grid_manager',
                        'dimension_type': '1d',
                        'dimension_name': 'n_levels'
                    })
                else:
                    # Surface emissions (common)
                    if config.has_size_bins or (hasattr(config, 'species') and len(config.species) > 1):
                        # Multi-species or size-resolved emissions
                        result.update({
                            'data_type': 'DIAG_REAL_3D',
                            'dimensions': ['nx', 'ny', 'n_species'],
                            'dimension_vars': ['dims_3d'],
                            'fortran_dims': 'dims_3d',
                            'dimension_source': 'process_config',
                            'dimension_type': '1d',
                            'dimension_name': 'n_species'
                        })
            elif config.process_type in ['chemistry', 'transport']:
                # Usually affect full column
                result.update({
                    'data_type': 'DIAG_REAL_3D',
                    'dimensions': ['nx', 'ny', 'nz'],
                    'dimension_vars': ['dims_3d_levels'],
                    'fortran_dims': 'dims_3d_levels',
                    'dimension_source': 'grid_manager',
                    'dimension_type': '1d',
                    'dimension_name': 'n_levels'
                })

        # 6. Infer from standard units
        elif '/m2/' in units:  # Surface flux units (per square meter)
            # Keep 2D inference (already set as default) - scalar per grid cell
            result.update({
                'dimension_type': 'scalar',
                'dimension_name': None
            })
        elif '/m3/' in units:  # Volume concentration units (per cubic meter)
            result.update({
                'data_type': 'DIAG_REAL_3D',
                'dimensions': ['nx', 'ny', 'nz'],
                'dimension_vars': ['dims_3d_levels'],
                'fortran_dims': 'dims_3d_levels',
                'dimension_source': 'grid_manager',
                'dimension_type': '1d',
                'dimension_name': 'n_levels'
            })

        # 7. Check diagnostic location (process vs scheme level)
        elif scheme_config:  # Scheme-specific diagnostic
            # Often more detailed (per-bin, per-level)
            if result['data_type'] == 'DIAG_REAL_2D':  # Upgrade to 3D if still 2D
                if config.has_size_bins or (hasattr(config, 'species') and len(config.species) > 1):
                    result.update({
                        'data_type': 'DIAG_REAL_3D',
                        'dimensions': ['nx', 'ny', 'n_species'],
                        'dimension_vars': ['dims_3d'],
                        'fortran_dims': 'dims_3d',
                        'dimension_source': 'process_config',
                        'dimension_type': '1d',
                        'dimension_name': 'n_species'
                    })

        return result

    def _analyze_required_dimensions(self, config: ProcessConfig) -> Dict[str, bool]:
        """Analyze which dimension arrays are actually needed for a process configuration."""
        required_dims = {
            'dims_2d': False,
            'dims_3d_species': False,
            'dims_3d_levels': False,
            'dims_soil': False,
            'dims_landuse': False,
            'dims_vegetation': False
        }

        # Always need 2D for most surface diagnostics
        required_dims['dims_2d'] = True

        # Check common diagnostics
        if hasattr(config, 'diagnostics') and config.diagnostics:
            for diagnostic in config.diagnostics:
                diag_props = self._infer_diagnostic_properties(diagnostic, config)
                fortran_dims = diag_props.get('fortran_dims', 'dims_2d')
                if fortran_dims in required_dims:
                    required_dims[fortran_dims] = True

        # Check scheme-specific diagnostics
        if hasattr(config, 'schemes') and config.schemes:
            for scheme in config.schemes:
                if hasattr(scheme, 'scheme_diagnostics') and scheme.scheme_diagnostics:
                    for diagnostic in scheme.scheme_diagnostics:
                        diag_props = self._infer_diagnostic_properties(diagnostic, config, scheme)
                        fortran_dims = diag_props.get('fortran_dims', 'dims_2d')
                        if fortran_dims in required_dims:
                            required_dims[fortran_dims] = True

        return required_dims

    def _get_all_met_fields_filter(self, scheme_dict: Dict[str, Any], config_dict: Dict[str, Any]) -> List[str]:
        """Jinja2 filter to get all required meteorological fields (common + scheme-specific)."""
        all_fields = set()

        # Add common process-level fields
        if 'required_met_fields' in config_dict and config_dict['required_met_fields']:
            all_fields.update(_extract_field_names(config_dict['required_met_fields']))

        # Add scheme-specific fields
        if 'required_met_fields' in scheme_dict and scheme_dict['required_met_fields']:
            all_fields.update(_extract_field_names(scheme_dict['required_met_fields']))

        # Return sorted list for consistent ordering
        return sorted(list(all_fields))

    def _format_fortran_dims(self, dimensions: List[str]) -> str:
        """Format dimension list for Fortran array declaration."""
        if len(dimensions) == 1:
            return 'dims_1d'
        elif len(dimensions) == 2:
            return 'dims_2d'
        elif len(dimensions) == 3:
            # Determine specific 3D dimension type based on third dimension
            third_dim = dimensions[2].lower()
            if 'nz' in third_dim or 'n_levels' in third_dim:
                return 'dims_3d_levels'
            elif 'n_soil' in third_dim:
                return 'dims_soil'
            elif 'n_landuse' in third_dim or 'n_veg' in third_dim:
                return 'dims_landuse'
            elif 'n_categories' in third_dim or 'n_types' in third_dim:
                return 'dims_3d'  # Default for species/bins
            else:
                return 'dims_3d'  # Default for species/bins
        elif len(dimensions) == 4:
            return 'dims_4d'
        else:
            return f"[{', '.join(dimensions)}]"

    def _get_dimension_source_and_access(self, dim_name: str) -> tuple:
        """
        Get the source manager and access pattern for a given dimension.

        Returns:
            tuple: (source_manager, access_pattern, variable_declaration)
        """

        # All dimensions come from GridManager for consistency
        grid_mapping = {
            'nx': ('grid_manager', 'nx', 'integer :: nx, ny, nz'),
            'ny': ('grid_manager', 'ny', 'integer :: nx, ny, nz'),
            'nz': ('grid_manager', 'nz', 'integer :: nx, ny, nz'),
            'nlev': ('grid_manager', 'nz', 'integer :: nx, ny, nz'),
        }

        if dim_name in grid_mapping:
            return grid_mapping[dim_name]
        else:
            # Default to process config for unknown dimensions
            return ('process_config', f'process_config%{dim_name}',
                   f'integer :: {dim_name}')

    def _generate_dimension_access_code(self, diagnostics: List[Dict], template_vars: Dict) -> str:
        """Generate code for accessing dimensional information."""

        required_dims = set()
        grid_needed = False

        # Collect all required dimensions
        for diag in diagnostics:
            if 'dimensions' in diag:
                for dim in diag['dimensions']:
                    if dim not in ['nx', 'ny']:  # Skip basic grid dims
                        required_dims.add(dim)
                        source, _, _ = self._get_dimension_source_and_access(dim)
                        if source == 'grid_manager':
                            grid_needed = True

        code_lines = []

        if grid_needed:
            code_lines.extend([
                "      ! Get grid manager pointer for dimension access",
                "      type(GridManagerType), pointer :: grid_mgr",
                "      grid_mgr => state_manager%get_grid_manager()",
                ""
            ])

            # Get basic grid dimensions
            code_lines.extend([
                "      ! Get basic grid dimensions",
                "      integer :: nx, ny, nz",
                "      call grid_mgr%get_dimensions(nx, ny, nz)",
                ""
            ])

        return '\n'.join(code_lines)

    def validate_config(self, config: ProcessConfig) -> None:
        """Validate process configuration.

        Args:
            config: Process configuration to validate

        Raises:
            ProcessValidationError: If configuration is invalid
        """
        errors = []

        # Basic validation
        if not config.name:
            errors.append("Process name is required")

        if not config.class_name:
            errors.append("Process class name is required")

        if not config.description:
            errors.append("Process description is required")

        if not config.author:
            errors.append("Author name is required")

        # Name validation
        if not config.name.replace('_', '').replace('-', '').isalnum():
            errors.append("Process name must be alphanumeric (with _ or - allowed)")

        if not config.class_name.replace('_', '').isalnum():
            errors.append("Class name must be alphanumeric (with _ allowed)")

        # Scheme validation
        if not config.schemes:
            errors.append("At least one scheme must be defined")

        scheme_names = [scheme.name for scheme in config.schemes]
        if len(set(scheme_names)) != len(scheme_names):
            errors.append("Scheme names must be unique")

        if config.default_scheme and config.default_scheme not in scheme_names:
            errors.append(f"Default scheme '{config.default_scheme}' not found in schemes")

        # Met field & species validation
        field_classifier = MetFieldClassification(self.metstate_file)
        known_met_fields = set(field_classifier.get_all_2d_fields() +
                               field_classifier.get_all_3d_atmospheric_fields() +
                               field_classifier.get_all_categorical_fields())

        all_met_fields = self.get_all_required_met_fields_combined(config)
        for field in all_met_fields:
            if known_met_fields and field not in known_met_fields:
                logger.warning(f"Unrecognized meteorological field '{field}' - ensure it is provided by MetState")

        if config.process_behavior and isinstance(config.process_behavior, dict):
            species_filter = config.process_behavior.get('species_filter', {})
            if species_filter.get('type') == 'by_metadata' and not species_filter.get('metadata_flags'):
                errors.append("species_filter type 'by_metadata' requires 'metadata_flags'")

        # Species validation
        if config.species and len(set(config.species)) != len(config.species):
            errors.append("Species names must be unique")

        # Output directory validation
        if config.output_dir:
            output_path = Path(config.output_dir)
            if output_path.exists() and not output_path.is_dir():
                errors.append(f"Output path exists but is not a directory: {config.output_dir}")

        if errors:
            raise ProcessValidationError("\n".join(errors))

    def has_persistent_state_variables(self, config: ProcessConfig) -> bool:
        """Check if any scheme has persistent state variables."""
        for scheme in config.schemes:
            if scheme.persistent_state_variables:
                return True
        return False

    def get_all_persistent_state_variables(self, config: ProcessConfig) -> Dict[str, List[Dict[str, Any]]]:
        """Get all persistent state variables organized by scheme."""
        all_variables = {}
        for scheme in config.schemes:
            if scheme.persistent_state_variables:
                processed_vars = []
                for var in scheme.persistent_state_variables:
                    processed_var = var.copy()
                    processed_var['clean_name'] = self._clean_state_variable_name(var['name'])
                    processed_var['fortran_type'] = self._infer_state_variable_type(var.get('default'), var['name'])
                    processed_var['dimensions'] = self._get_state_variable_dimensions(var['name'])
                    processed_var['is_allocatable'] = '(:)' in var['name']
                    processed_vars.append(processed_var)
                all_variables[scheme.name] = processed_vars
        return all_variables

    def load_config(self, config_path: Union[str, Path]) -> ProcessConfig:
        """Load and validate process configuration from YAML file.

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Validated ProcessConfig object

        Raises:
            ProcessValidationError: If configuration is invalid
            FileNotFoundError: If config file doesn't exist
            yaml.YAMLError: If YAML parsing fails
        """
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, 'r') as f:
            data = yaml.safe_load(f)

        def _make_dataclass(cls, kwargs_dict):
            if not isinstance(kwargs_dict, dict):
                return kwargs_dict
            valid_keys = {field.name for field in dataclasses.fields(cls)}
            return cls(**{k: v for k, v in kwargs_dict.items() if k in valid_keys})

        # Convert process_behavior to ProcessBehavior object if present
        if 'process_behavior' in data and data['process_behavior']:
            data['process_behavior'] = _make_dataclass(ProcessBehavior, data['process_behavior'])

        # Convert schemes to SchemeConfig objects
        schemes = []
        schemes_data = data.get('schemes', {})
        if isinstance(schemes_data, dict):
            # Schemes defined as dict with keys
            for scheme_name, scheme_data in schemes_data.items():
                # Add the scheme name if not already present
                if 'name' not in scheme_data:
                    scheme_data['name'] = scheme_name

                # Convert scheme_behavior to SchemeBehavior object if present
                if 'scheme_behavior' in scheme_data and scheme_data['scheme_behavior']:
                    scheme_data['scheme_behavior'] = _make_dataclass(SchemeBehavior, scheme_data['scheme_behavior'])

                # Handle diagnostics field name change
                if 'diagnostics' in scheme_data:
                    scheme_data['scheme_diagnostics'] = scheme_data.pop('diagnostics')

                # Add gas_or_aero field if present, default to 'both'
                if 'gas_or_aero' not in scheme_data:
                    scheme_data['gas_or_aero'] = 'both'

                scheme = _make_dataclass(SchemeConfig, scheme_data)
                schemes.append(scheme)
        elif isinstance(schemes_data, list):
            # Schemes defined as list
            for scheme_data in schemes_data:
                # Convert scheme_behavior to SchemeBehavior object if present
                if 'scheme_behavior' in scheme_data and scheme_data['scheme_behavior']:
                    scheme_data['scheme_behavior'] = _make_dataclass(SchemeBehavior, scheme_data['scheme_behavior'])

                # Handle diagnostics field name change
                if 'diagnostics' in scheme_data:
                    scheme_data['scheme_diagnostics'] = scheme_data.pop('diagnostics')

                # Add gas_or_aero field if present, default to 'both'
                if 'gas_or_aero' not in scheme_data:
                    scheme_data['gas_or_aero'] = 'both'

                scheme = _make_dataclass(SchemeConfig, scheme_data)
                schemes.append(scheme)

        # Remove schemes from data and create ProcessConfig
        data['schemes'] = schemes
        config = _make_dataclass(ProcessConfig, data)

        # Load species database if process behavior specifies species filtering
        if config.process_behavior and config.process_behavior.species_filter:
            config.species = self._load_filtered_species(config.process_behavior.species_filter)

        # Validate configuration
        self.validate_config(config)

        return config

    def get_all_required_species_properties(self, config: ProcessConfig) -> List[str]:
        """Collect all unique required species properties from all schemes.

        Args:
            config: ProcessConfig object containing schemes

        Returns:
            List of unique property names required by all schemes
        """
        all_properties = set()
        for scheme in config.schemes:
            if scheme.required_species_properties:
                all_properties.update(scheme.required_species_properties)

        # Return sorted list for consistent ordering
        return sorted(list(all_properties))

    def get_all_required_met_fields(self, config: ProcessConfig, scheme_config: SchemeConfig = None) -> List[str]:
        """Get all unique required meteorological fields for a scheme (common + scheme-specific).

        Args:
            config: ProcessConfig object containing common fields
            scheme_config: Optional SchemeConfig object containing scheme-specific fields

        Returns:
            List of unique field names required by the process and/or scheme
        """
        all_fields = set()

        # Add common process-level fields
        if hasattr(config, 'required_met_fields') and config.required_met_fields:
            all_fields.update(_extract_field_names(config.required_met_fields))

        # Add scheme-specific fields
        if scheme_config and hasattr(scheme_config, 'required_met_fields') and scheme_config.required_met_fields:
            all_fields.update(_extract_field_names(scheme_config.required_met_fields))

        # Return sorted list for consistent ordering
        return sorted(list(all_fields))

    def get_gas_schemes(self, config: ProcessConfig) -> List[SchemeConfig]:
        """Get schemes that apply to gas species."""
        return [scheme for scheme in config.schemes if scheme.gas_or_aero in ['gas', 'both']]

    def get_aero_schemes(self, config: ProcessConfig) -> List[SchemeConfig]:
        """Get schemes that apply to aerosol species."""
        return [scheme for scheme in config.schemes if scheme.gas_or_aero in ['aero', 'both']]

    def get_all_required_constants(self, config: ProcessConfig) -> List[str]:
        """Collect all unique required constants from all schemes.

        Args:
            config: ProcessConfig object containing schemes

        Returns:
            List of unique constant names required by all schemes
        """
        all_constants = set()
        # Collect process-level required constants
        if config.required_constants:
            all_constants.update(config.required_constants)
        # Collect scheme-level required constants
        for scheme in config.schemes:
            if scheme.required_constants:
                all_constants.update(scheme.required_constants)

        # Return sorted list for consistent ordering
        return sorted(list(all_constants))

    def has_required_time_parameters(self, config: ProcessConfig) -> bool:
        """Check if any scheme requires time parameters.

        Args:
            config: ProcessConfig object containing schemes

        Returns:
            True if any scheme has required_time_parameters, False otherwise
        """
        for scheme in config.schemes:
            if scheme.required_time_parameters:
                return True
        return False

    def get_all_required_time_parameters(self, config: ProcessConfig) -> List[str]:
        """Collect all unique required time parameters from all schemes.

        Args:
            config: ProcessConfig object containing schemes

        Returns:
            List of unique time parameter names required by all schemes
        """
        all_time_params = set()
        for scheme in config.schemes:
            if scheme.required_time_parameters:
                all_time_params.update(scheme.required_time_parameters)

        # Return sorted list for consistent ordering
        return sorted(list(all_time_params))

    def has_persistent_state_variables(self, config: ProcessConfig) -> bool:
        """Check if any scheme has persistent state variables."""
        for scheme in config.schemes:
            if scheme.persistent_state_variables:
                return True
        return False

    def get_all_persistent_state_variables(self, config: ProcessConfig) -> Dict[str, List[Dict[str, Any]]]:
        """Get all persistent state variables organized by scheme."""
        all_variables = {}
        for scheme in config.schemes:
            if scheme.persistent_state_variables:
                processed_vars = []
                for var in scheme.persistent_state_variables:
                    processed_var = var.copy()
                    processed_var['clean_name'] = self._clean_state_variable_name(var['name'])
                    processed_var['fortran_type'] = self._infer_state_variable_type(var.get('default'), var['name'])
                    processed_var['dimensions'] = self._get_state_variable_dimensions(var['name'])
                    processed_var['is_allocatable'] = '(:)' in var['name']
                    processed_vars.append(processed_var)
                all_variables[scheme.name] = processed_vars
        return all_variables

    def _load_filtered_species(self, species_filter: Dict[str, Any]) -> List[str]:
        """Load species based on filter criteria.

        Args:
            species_filter: Dictionary specifying filter criteria

        Returns:
            List of species names that match the filter
        """
        filter_type = species_filter.get('type', 'all_species')

        if filter_type == 'all_species':
            return []  # Will be filled by ChemState

        elif filter_type == 'by_list':
            return species_filter.get('species_list', [])

        elif filter_type == 'by_metadata':
            # Try to load species database
            species_db_path = Path(self.template_dir.parent / "configs" / "species_database.yaml")
            if not species_db_path.exists():
                logger.warning(f"Species database not found at {species_db_path}, using empty species list")
                return []

            try:
                with open(species_db_path, 'r') as f:
                    species_db = yaml.safe_load(f)

                metadata_flags = species_filter.get('metadata_flags', [])
                filtered_species = []

                for species_name, species_data in species_db.get('species_database', {}).items():
                    # Check if species has all required metadata flags set to True
                    if all(species_data.get(flag, False) for flag in metadata_flags):
                        filtered_species.append(species_name)

                logger.info(f"Filtered {len(filtered_species)} species using metadata flags: {metadata_flags}")
                return filtered_species

            except Exception as e:
                logger.warning(f"Error loading species database: {e}, using empty species list")
                return []

        elif filter_type == 'emission_mapping':
            # Load species from emission mapping configuration
            # The species will be loaded dynamically at runtime from ConfigManager emission mapping
            # Return empty list here - actual species loading handled by load_species_from_emission_mapping
            logger.info("Using emission_mapping species filter - species will be loaded from ConfigManager at runtime")
            return []

        return []

    def get_all_required_met_fields_combined(self, config: ProcessConfig) -> List[str]:
        """Collect all unique required met fields across process and all schemes."""
        all_fields = set()
        if hasattr(config, 'required_met_fields') and config.required_met_fields:
            all_fields.update(_extract_field_names(config.required_met_fields))
        for scheme in config.schemes:
            if hasattr(scheme, 'required_met_fields') and scheme.required_met_fields:
                all_fields.update(_extract_field_names(scheme.required_met_fields))
            elif isinstance(scheme, dict) and 'required_met_fields' in scheme:
                all_fields.update(_extract_field_names(scheme['required_met_fields']))
        return sorted(list(all_fields))

    def generate_process(self, config: ProcessConfig) -> None:
        """Generate complete process implementation.

        Args:
            config: Validated process configuration
        """
        logger.info(f"Generating process: {config.name}")

        script_dir = Path(__file__).resolve().parent
        repo_root = script_dir.parent.parent

        if config.output_dir:
            process_dir = Path(config.output_dir)
            test_dir = process_dir / "tests"
            docs_dir = process_dir
        else:
            process_dir = repo_root / config.src_base_dir / config.name
            test_dir = repo_root / "tests" / "process" / config.name
            docs_dir = repo_root / "docs" / "processes" / config.name

        # Create directory structure
        self._create_directory_structure(process_dir, test_dir, docs_dir, config)

        # Generate files
        self._generate_cpp_wrapper(process_dir, config)
        self._generate_science_bridge(process_dir, config)
        self._generate_common_module(process_dir, config)
        self._generate_schemes(process_dir, config)
        self._generate_cmake_files(process_dir, config)

        if config.generate_tests:
            self._generate_tests(test_dir, config)

        if config.generate_docs:
            self._generate_documentation(docs_dir, config)

        logger.info(f"Process generation complete: {process_dir}")

    def _create_directory_structure(self, process_dir: Path, test_dir: Path, docs_dir: Path, config: ProcessConfig) -> None:
        """Create the directory structure for the process."""
        logger.info(f"Creating directory structure: {process_dir}")

        directories = [
            process_dir,
            process_dir / "schemes",
        ]

        if config.generate_tests:
            directories.append(test_dir)

        if config.generate_docs and docs_dir != process_dir:
            directories.append(docs_dir)

        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)

    def _generate_cpp_wrapper(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate C++ Process Wrapper header and source files."""
        logger.info("Generating C++ process wrapper header and source")

        all_met_fields = self.get_all_required_met_fields_combined(config)

        # 1. Header
        template_hpp = self.env.get_template('catchem_process.hpp.j2')
        content_hpp = template_hpp.render(
            config=config,
            all_met_fields=all_met_fields,
            timestamp=datetime.now().isoformat()
        )
        file_hpp = process_dir / f"catchem_process_{config.name}.hpp"
        with open(file_hpp, 'w') as f:
            f.write(content_hpp)
        logger.info(f"Generated: {file_hpp}")

        # 2. Source
        template_cpp = self.env.get_template('catchem_process.cpp.j2')
        content_cpp = template_cpp.render(
            config=config,
            all_met_fields=all_met_fields,
            timestamp=datetime.now().isoformat()
        )
        file_cpp = process_dir / f"catchem_process_{config.name}.cpp"
        with open(file_cpp, 'w') as f:
            f.write(content_cpp)
        logger.info(f"Generated: {file_cpp}")

    def _generate_science_bridge(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate Fortran BIND(C) Science Bridge module."""
        logger.info("Generating Fortran science bridge")

        all_met_fields = self.get_all_required_met_fields_combined(config)

        template = self.env.get_template('science_bridge.F90.j2')
        content = template.render(
            config=config,
            all_met_fields=all_met_fields,
            timestamp=datetime.now().isoformat()
        )

        filename = f"{config.class_name}ScienceBridge.F90"
        output_file = process_dir / filename

        with open(output_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {output_file}")

    def _generate_common_module(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate the common types and utilities module."""
        logger.info("Generating common module")

        template = self.env.get_template('process_common.F90.j2')

        all_required_species_properties = self.get_all_required_species_properties(config)
        field_classifier = MetFieldClassification(self.metstate_file)

        content = template.render(
            config=config,
            all_required_species_properties=all_required_species_properties,
            field_classifier=field_classifier,
            has_persistent_state_variables=self.has_persistent_state_variables(config),
            all_persistent_state_variables=self.get_all_persistent_state_variables(config),
            generation_date=datetime.now().isoformat(),
            version=config.version,
            timestamp=datetime.now().isoformat(),
            gas_schemes=self.get_gas_schemes(config),
            aero_schemes=self.get_aero_schemes(config)
        )

        filename = f"{config.class_name}Common_Mod.F90"
        output_file = process_dir / filename

        with open(output_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {output_file}")

    def _generate_schemes(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate scheme implementation modules."""
        logger.info("Generating scheme modules")

        from dataclasses import asdict, is_dataclass
        def to_dict(obj):
            if isinstance(obj, dict):
                return obj
            elif is_dataclass(obj):
                return asdict(obj)
            else:
                return obj.__dict__

        field_classifier = MetFieldClassification(self.metstate_file)
        template = self.env.get_template('scheme_module.F90.j2')
        schemes_dir = process_dir / "schemes"

        for scheme in config.schemes:
            scheme_dict = to_dict(scheme)
            config_dict = to_dict(config)

            self.env.globals['config'] = config

            try:
                content = template.render(
                    config=config_dict,
                    scheme=scheme_dict,
                    all_required_constants=self.get_all_required_constants(config),
                    needs_time_state=self.has_required_time_parameters(config),
                    all_required_time_parameters=self.get_all_required_time_parameters(config),
                    field_classifier=field_classifier,
                    has_persistent_state_variables=self.has_persistent_state_variables(config),
                    all_persistent_state_variables=self.get_all_persistent_state_variables(config),
                    timestamp=datetime.now().isoformat()
                )
            except Exception as e:
                logger.error(f"Template rendering failed for scheme {scheme_dict.get('name', '<unknown>')}: {e}")
                content = f"! Template rendering failed: {e}\n"

            filename = f"{config_dict['class_name']}Scheme_{scheme_dict['class_name']}_Mod.F90"
            output_file = schemes_dir / filename

            with open(output_file, 'w') as f:
                f.write(content)

            logger.info(f"Generated scheme: {output_file}")

    def _generate_cmake_files(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate CMake configuration files."""
        logger.info("Generating CMake files")

        template = self.env.get_template('CMakeLists.txt.j2')
        content = template.render(config=config, timestamp=datetime.now().isoformat())

        cmake_file = process_dir / "CMakeLists.txt"
        with open(cmake_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {cmake_file}")

        schemes_cmake_template = self.env.get_template('schemes_CMakeLists.txt.j2')
        schemes_content = schemes_cmake_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        schemes_cmake_file = process_dir / "schemes" / "CMakeLists.txt"
        with open(schemes_cmake_file, 'w') as f:
            f.write(schemes_content)

        logger.info(f"Generated: {schemes_cmake_file}")

    def _generate_tests(self, test_dir: Path, config: ProcessConfig) -> None:
        """Generate standalone Fortran CTest science test."""
        logger.info(f"Generating test files in: {test_dir}")

        all_met_fields = self.get_all_required_met_fields_combined(config)

        test_dir.mkdir(parents=True, exist_ok=True)

        science_test_template = self.env.get_template('test_science.f90.j2')
        content = science_test_template.render(
            config=config,
            all_met_fields=all_met_fields,
            timestamp=datetime.now().isoformat()
        )

        test_file = test_dir / f"test_{config.name}_science.f90"
        with open(test_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated science test: {test_file}")

        test_cmake_template = self.env.get_template('test_CMakeLists.txt.j2')
        test_cmake_content = test_cmake_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        test_dir.mkdir(parents=True, exist_ok=True)
        test_cmake_file = test_dir / "CMakeLists.txt"
        with open(test_cmake_file, 'w') as f:
            f.write(test_cmake_content)

    def _generate_documentation(self, docs_dir: Path, config: ProcessConfig) -> None:
        """Generate consolidated documentation in docs/processes/<process_name> directory."""
        logger.info(f"Generating documentation in: {docs_dir}")

        # Single comprehensive documentation file
        doc_template = self.env.get_template('process_documentation.md.j2')
        doc_content = doc_template.render(config=config, timestamp=datetime.now().isoformat())

        doc_file = docs_dir / f"{config.name}.md"
        with open(doc_file, 'w') as f:
            f.write(doc_content)

        logger.info(f"Generated documentation in: {docs_dir}")



    def generate_template_config(self, process_type: str = "emission") -> Dict[str, Any]:
        """Generate a template configuration for a given process type.

        Args:
            process_type: Type of process (emission, chemistry, transport, etc.)

        Returns:
            Dictionary containing template configuration
        """
        template_configs = {
            "emission": {
                "name": "my_emission",
                "class_name": "MyEmission",
                "description": "Description of my emission process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "emission",
                "process_behavior": {
                    "type": "source",
                    "tendency_mode": "additive",
                    "parallelization": "column",
                    "spatial_scope": "column",
                    "timestep_dependency": "independent"
                },
                "is_multiphase": False,
                "has_size_bins": False,
                "species": ["species1", "species2"],
                "schemes": [
                    {
                        "name": "simple",
                        "class_name": "Simple",
                        "description": "Simple emission scheme",
                        "author": "Your Name",
                        "required_met_fields": ["temperature", "wind_speed"],
                        "diagnostics": [
                            {
                                "name": "emission_rate",
                                "units": "kg/m2/s",
                                "description": "Emission rate"
                            }
                        ]
                    }
                ],
                "default_scheme": "simple",
                "required_met_fields": ["temperature", "wind_speed"],
                "diagnostics": [
                    {
                        "name": "total_emissions",
                        "units": "kg/m2/s",
                        "description": "Total emission rate"
                    }
                ],
                "generate_tests": True,
                "generate_docs": True
            },

            "chemistry": {
                "name": "my_chemistry",
                "class_name": "MyChemistry",
                "description": "Description of my chemistry process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "chemistry",
                "process_behavior": {
                    "type": "transformation",
                    "tendency_mode": "replacement",
                    "parallelization": "column",
                    "spatial_scope": "column",
                    "timestep_dependency": "dependent"
                },
                "is_multiphase": True,
                "has_size_bins": False,
                "species": ["O3", "NO", "NO2", "OH"],
                "phases": ["gas", "aqueous"],
                "schemes": [
                    {
                        "name": "mechanism1",
                        "class_name": "Mechanism1",
                        "description": "Chemistry mechanism 1",
                        "author": "Your Name",
                        "algorithm_type": "implicit",
                        "required_met_fields": ["temperature", "pressure", "humidity"],
                        "diagnostics": [
                            {
                                "name": "reaction_rates",
                                "units": "molec/cm3/s",
                                "description": "Chemical reaction rates"
                            }
                        ]
                    }
                ],
                "default_scheme": "mechanism1",
                "required_met_fields": ["temperature", "pressure", "humidity"],
                "timestep_dependency": "dependent",
                "parallelization": "column",
                "generate_tests": True,
                "generate_docs": True
            },

            "deposition": {
                "name": "my_deposition",
                "class_name": "MyDeposition",
                "description": "Dry and wet deposition process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "deposition",
                "process_behavior": {
                    "type": "sink",
                    "tendency_mode": "additive",
                    "parallelization": "column",
                    "spatial_scope": "column"
                },
                "species": ["O3", "NO2", "SO2"],
                "schemes": [
                    {
                        "name": "wesely",
                        "class_name": "Wesely",
                        "description": "Wesely dry deposition scheme",
                        "author": "Your Name",
                        "required_met_fields": ["USTAR", "Z0", "TSKIN"]
                    }
                ],
                "default_scheme": "wesely",
                "required_met_fields": ["USTAR", "Z0", "TSKIN"],
                "generate_tests": True,
                "generate_docs": True
            },

            "settling": {
                "name": "my_settling",
                "class_name": "MySettling",
                "description": "Gravitational settling process for aerosols",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "settling",
                "process_behavior": {
                    "type": "sink",
                    "tendency_mode": "additive",
                    "parallelization": "column",
                    "spatial_scope": "column"
                },
                "species": ["DUST1", "DUST2", "SEASALT1"],
                "schemes": [
                    {
                        "name": "stokes",
                        "class_name": "Stokes",
                        "description": "Stokes gravitational settling scheme",
                        "author": "Your Name",
                        "required_met_fields": ["AIRDEN", "T", "P"]
                    }
                ],
                "default_scheme": "stokes",
                "required_met_fields": ["AIRDEN", "T", "P"],
                "generate_tests": True,
                "generate_docs": True
            },

            "plume_rise": {
                "name": "my_plume_rise",
                "class_name": "MyPlumeRise",
                "description": "Wildfire and point source plume rise process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "plume_rise",
                "process_behavior": {
                    "type": "transport",
                    "tendency_mode": "additive",
                    "parallelization": "column",
                    "spatial_scope": "column"
                },
                "species": ["CO", "NO2", "PM25"],
                "schemes": [
                    {
                        "name": "freitas",
                        "class_name": "Freitas",
                        "description": "Freitas 1D plume rise scheme",
                        "author": "Your Name",
                        "required_met_fields": ["T", "PBLH", "U10M", "V10M"]
                    }
                ],
                "default_scheme": "freitas",
                "required_met_fields": ["T", "PBLH", "U10M", "V10M"],
                "generate_tests": True,
                "generate_docs": True
            }
        }

        return template_configs.get(process_type, template_configs["emission"])


def main():
    """Main entry point for the process generator."""
    parser = argparse.ArgumentParser(
        description="CATChem Process Generator Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate a process from configuration
  python process_generator.py generate --config my_process.yaml

  # Generate with automatic MetState field discovery
  python process_generator.py generate --config my_process.yaml --metstate src/core/metstate_mod.F90

  # Validate a configuration file
  python process_generator.py validate --config my_process.yaml

  # Generate template configuration
  python process_generator.py template --type emission --output emission_template.yaml

  # Generate with custom template directory
  python process_generator.py generate --config my_process.yaml --templates ./my_templates

  # Inspect discovered meteorological fields
  python process_generator.py fields --type all --verbose

  # Show specific field types
  python process_generator.py fields --type 2d --verbose
  python process_generator.py fields --type 3d --verbose

Features:
  # Automatic MetState Field Discovery
  The generator automatically discovers meteorological fields from the MetState definition file.
  This ensures that field classifications (2D surface, 3D atmospheric, categorical) are always
  up-to-date with the actual MetState implementation. The generator will attempt to find the
  MetState file automatically, or you can specify it explicitly with --metstate.

  Field Types:
  - 2D Surface: FROCEAN, SST, USTAR, PBLH (scalar access)
  - 3D Atmospheric: T, QV, P, U, V (array access when affects_full_column=true)
  - Categorical: SOILM, FRLANDUSE (special dimension arrays)
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Generate command
    generate_parser = subparsers.add_parser('generate',
                                           help='Generate process implementation',
                                           description='Generate complete process implementation from YAML configuration. '
                                                     'Renders modern C++ process classes and Fortran science bridges.')
    generate_parser.add_argument('--config', '-c', required=True,
                               help='Path to YAML configuration file')
    generate_parser.add_argument('--output-dir', '--output', '-o', dest='output',
                               help='Output directory for generated files (default: src/process/<process_name>)')
    generate_parser.add_argument('--templates', '-t',
                               help='Path to template directory')
    generate_parser.add_argument('--metstate-file', '--metstate', '-m', dest='metstate',
                               help='Path to MetState file for automatic field discovery')
    generate_parser.add_argument('--force', '-f', action='store_true',
                               help='Force overwriting existing target files without prompting')
    generate_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose output')

    # Validate command
    validate_parser = subparsers.add_parser('validate',
                                           help='Validate configuration file',
                                           description='Validate YAML configuration file syntax and check field compatibility. '
                                                     'Ensures configuration is ready for generation.')
    validate_parser.add_argument('--config', '-c', required=True,
                               help='Path to YAML configuration file')
    validate_parser.add_argument('--metstate-file', '--metstate', '-m', dest='metstate',
                               help='Path to MetState file for automatic field discovery')
    validate_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose output')

    # Template command
    template_parser = subparsers.add_parser('template', help='Generate template configuration')
    template_parser.add_argument('--type', '-t', choices=['emission', 'chemistry', 'deposition', 'settling', 'plume_rise'],
                               default='emission', help='Type of process template')
    template_parser.add_argument('--output', '-o', required=True,
                               help='Output file for template')

    # Fields command - show discovered MetState fields
    fields_parser = subparsers.add_parser('fields',
                                         help='Show discovered MetState fields',
                                         description='Display meteorological fields discovered from MetState definition. '
                                                   'Fields are automatically classified into 2D surface, 3D atmospheric, '
                                                   'and categorical types for proper code generation.')
    fields_parser.add_argument('--metstate', '-m',
                               help='Path to MetState file for field discovery')
    fields_parser.add_argument('--type', '-t', choices=['2d', '3d', 'categorical', 'all'],
                               default='all', help='Type of fields to show')
    fields_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose output')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Set up logging level
    if hasattr(args, 'verbose') and args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        if args.command == 'generate':
            generator = ProcessGenerator(args.templates, getattr(args, 'metstate', None))
            config = generator.load_config(args.config)

            # Override output directory if specified on command line
            if hasattr(args, 'output') and args.output:
                config.output_dir = args.output

            generator.generate_process(config)

        elif args.command == 'validate':
            generator = ProcessGenerator(metstate_file=getattr(args, 'metstate', None))
            config = generator.load_config(args.config)
            print(f"Configuration '{Path(args.config).name}' is valid.")
            logger.info(f"Configuration is valid: {args.config}")
            return 0

        elif args.command == 'template':
            generator = ProcessGenerator()
            template_config = generator.generate_template_config(args.type)

            with open(args.output, 'w') as f:
                yaml.dump(template_config, f, default_flow_style=False, sort_keys=False)

            logger.info(f"Template configuration written to: {args.output}")

        elif args.command == 'fields':
            generator = ProcessGenerator(metstate_file=getattr(args, 'metstate', None))
            classifier = MetFieldClassification(generator.metstate_file)

            if args.type in ['2d', 'all']:
                fields_2d = classifier.get_all_2d_fields()
                logger.info(f"2D Surface fields ({len(fields_2d)}): {', '.join(fields_2d[:10])}{'...' if len(fields_2d) > 10 else ''}")

            if args.type in ['3d', 'all']:
                fields_3d = classifier.get_all_3d_atmospheric_fields()
                logger.info(f"3D Atmospheric fields ({len(fields_3d)}): {', '.join(fields_3d[:10])}{'...' if len(fields_3d) > 10 else ''}")

            if args.type in ['categorical', 'all']:
                fields_cat = classifier.get_all_categorical_fields()
                logger.info(f"Categorical fields ({len(fields_cat)}): {', '.join(fields_cat)}")

            if args.verbose:
                if args.type in ['2d', 'all']:
                    print(f"\nAll 2D Surface fields ({len(fields_2d)}):")
                    for i, field in enumerate(fields_2d, 1):
                        print(f"  {i:3d}. {field}")

                if args.type in ['3d', 'all']:
                    print(f"\nAll 3D Atmospheric fields ({len(fields_3d)}):")
                    for i, field in enumerate(fields_3d, 1):
                        print(f"  {i:3d}. {field}")

                if args.type in ['categorical', 'all']:
                    print(f"\nAll Categorical fields ({len(fields_cat)}):")
                    for i, field in enumerate(fields_cat, 1):
                        print(f"  {i:3d}. {field}")

    except ProcessValidationError as e:
        logger.error(f"Configuration validation failed: {e}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 2
    except (IOError, OSError) as e:
        logger.error(f"I/O error: {e}")
        return 2
    except Exception as e:
        logger.error(f"Error: {e}")
        if hasattr(args, 'verbose') and args.verbose:
            import traceback
            traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
