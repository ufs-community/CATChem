#!/usr/bin/env python3
"""
Generate Fortran macros for MetStateType fields.

This script parses a Fortran source file containing the MetStateType definition and generates Fortran include files
for field accessors, allocation, and deallocation macros. The generated macros are used to implement field-by-field
allocation and access in a modular, robust way.

Usage
-----
Run from the command line:
    generate_metstate_macros.py accessor <input_f90> <output_inc>
    generate_metstate_macros.py allocate <input_f90> <output_inc>
    generate_metstate_macros.py deallocate <input_f90> <output_inc>
    generate_metstate_macros.py column_accessor <input_f90> <output_inc>
    generate_metstate_macros.py 2d_scalar_accessor <input_f90> <output_inc>
    generate_metstate_macros.py scalar_accessor <input_f90> <output_inc>
    generate_metstate_macros.py virtualcolumn_populate <input_f90> <output_inc>
    generate_metstate_macros.py virtualmet_populate <input_f90> <output_inc>

Parameters
----------
mode : str
    Which macro to generate. One of: accessor, allocate, deallocate, column_accessor, 2d_scalar_accessor, scalar_accessor, virtualcolumn_populate, virtualmet_populate.
input_f90 : str
    Path to the Fortran source file containing the MetStateType definition.
output_inc : str
    Path to the output .inc file to write.

Examples
--------
$ python generate_metstate_macros.py allocate ../src/core/metstate_mod.F90 metstate_allocate_fields.inc

"""
import re
import sys

def parse_metstate_type(filename):
    """
    Parse the MetStateType definition in a Fortran source file.

    Parameters
    ----------
    filename : str
        Path to the Fortran source file containing MetStateType.

    Returns
    -------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each allocatable field.
        type_name is 'real', 'integer', 'logical', or 'character'
        is_edge is True if the field uses nz+1 dimension (edge fields).
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    in_type = False
    fields = []  # List of (name, type_name, rank, dims, is_edge)
    for line in lines:
        if 'TYPE, PUBLIC :: MetStateType' in line:
            in_type = True
        elif in_type and 'end type' in line.lower():
            break
        elif in_type:
            # Match different field types: real(fp), integer, logical, character
            # REAL(fp), ALLOCATABLE :: name(dimensions) or name
            m_real_alloc = re.match(r'\s*REAL\(fp\),\s*ALLOCATABLE\s*::\s*(\w+)(?:\s*\(([^)]*)\))?\s*', line, re.IGNORECASE)
            # REAL(fp) :: name (scalar)
            m_real_scalar = re.match(r'\s*REAL\(fp\)\s*::\s*(\w+)(?:\s*=.*?)?\s*', line, re.IGNORECASE)
            # INTEGER, ALLOCATABLE :: name(dimensions) or name
            m_int_alloc = re.match(r'\s*INTEGER,\s*ALLOCATABLE\s*::\s*(\w+)(?:\s*\(([^)]*)\))?\s*', line, re.IGNORECASE)
            # INTEGER :: name (scalar)
            m_int_scalar = re.match(r'\s*INTEGER\s*::\s*(\w+)(?:\s*=.*?)?\s*', line, re.IGNORECASE)
            # LOGICAL, ALLOCATABLE :: name(dimensions) or name
            m_log_alloc = re.match(r'\s*LOGICAL,\s*ALLOCATABLE\s*::\s*(\w+)(?:\s*\(([^)]*)\))?\s*', line, re.IGNORECASE)
            # LOGICAL :: name (scalar)
            m_log_scalar = re.match(r'\s*LOGICAL\s*::\s*(\w+)(?:\s*=.*?)?\s*', line, re.IGNORECASE)
            # CHARACTER(len=*), ALLOCATABLE :: name(dimensions) or name
            m_char_alloc = re.match(r'\s*CHARACTER\s*\([^)]*\),\s*ALLOCATABLE\s*::\s*(\w+)(?:\s*\(([^)]*)\))?\s*', line, re.IGNORECASE)
            # CHARACTER(len=*) :: name (scalar)
            m_char_scalar = re.match(r'\s*CHARACTER\s*\([^)]*\)\s*::\s*(\w+)(?:\s*=.*?)?\s*', line, re.IGNORECASE)

            match = None
            type_name = None
            dims = None

            if m_real_alloc:
                match = m_real_alloc
                type_name = 'real'
                dims = match.group(2) if match.group(2) else None
            elif m_real_scalar:
                match = m_real_scalar
                type_name = 'real'
                dims = None
            elif m_int_alloc:
                match = m_int_alloc
                type_name = 'integer'
                dims = match.group(2) if match.group(2) else None
            elif m_int_scalar:
                match = m_int_scalar
                type_name = 'integer'
                dims = None
            elif m_log_alloc:
                match = m_log_alloc
                type_name = 'logical'
                dims = match.group(2) if match.group(2) else None
            elif m_log_scalar:
                match = m_log_scalar
                type_name = 'logical'
                dims = None
            elif m_char_alloc:
                match = m_char_alloc
                type_name = 'character'
                dims = match.group(2) if match.group(2) else None
            elif m_char_scalar:
                match = m_char_scalar
                type_name = 'character'
                dims = None

            if match:
                name = match.group(1)
                if dims:
                    rank = dims.count(',') + 1
                else:
                    rank = 0  # Scalar field
                # Check if this is an edge field (nz+1 dimension)
                is_edge = 'nz+1' in line or 'nlevs+1' in line.lower()
                fields.append((name, type_name, rank, dims, is_edge))

    return fields

def write_accessor(fields, output_file):
    """
    Write a macro for field accessors (all 2D/3D fields).

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Deduplicate by case-insensitive field name
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        key = name.lower()
        if key not in unique_fields:
            unique_fields[key] = (name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case(" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) then\n")
            if rank == 3:
                f.write(f"      column_ptr => this%{name}(col_i, col_j, :)\n")
            elif rank == 2:
                f.write(f"      column_ptr => null()\n")
            f.write(f"   endif\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_allocate(fields, output_file):
    """
    Write a macro for allocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.
    Handles edge fields (nz+1) properly.
    Now uses automatic detection of conditional allocation patterns.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 0:
                continue

            # Check if this field has conditional allocation
            conditional_info = get_conditional_allocation_info(name)

            if conditional_info:
                # Handle conditional allocation based on the field type
                if conditional_info['type'] == 'soil':
                    if rank == 3:
                        f.write(f"  if (nsoil > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif conditional_info['type'] == 'soiltype':
                    if rank == 3:
                        f.write(f"  if (nsoiltype > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif conditional_info['type'] == 'surftype':
                    if rank == 3:
                        f.write(f"  if (nSURFTYPE > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                # Standard allocation for non-conditional fields
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    if is_edge:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz+1))  ! Edge field\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")

        # Individual field cases
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")

            # Check if this field has conditional allocation
            conditional_info = get_conditional_allocation_info(name)

            if conditional_info:
                # Handle conditional allocation based on the field type
                if conditional_info['type'] == 'soil':
                    if rank == 3:
                        f.write(f"  if (nsoil > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif conditional_info['type'] == 'soiltype':
                    if rank == 3:
                        f.write(f"  if (nsoiltype > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif conditional_info['type'] == 'surftype':
                    if rank == 3:
                        f.write(f"  if (nSURFTYPE > 0) then\n")
                        f.write(f"    if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                        f.write(f"  endif\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                # Standard allocation for non-conditional fields
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    if is_edge:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz+1))  ! Edge field\n")
                    else:
                        f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")

        f.write("case default\n")
        f.write("  ! No allocation for unknown field\n")

def write_deallocate(fields, output_file):
    """
    Write a macro for deallocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 0:
                continue
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        # Individual field cases
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        f.write("case default\n")
        f.write("  ! No deallocation for unknown field\n")

def write_column_accessor(fields, output_file):
    """
    Write a macro for 3D column field accessors.
    Only includes REAL type fields since the interface returns real(fp) pointers.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Only 3D REAL fields
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 3 and type_name.upper() == 'REAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) column_ptr => this%{name}(col_i, col_j, :)\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def classify_fields(fields):
    """
    Classify MetState fields by type and suitability for VirtualMet.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.

    Returns
    -------
    tuple of lists
        (atmospheric_3d, categorical_3d, surface_2d, scalar_0d)
        where each list contains field names suitable for VirtualMet
    """
    # Categorical 3D fields with non-vertical dimensions (soil, land use, etc.)
    categorical_patterns = {
        'SOILM', 'SOILT', 'FRSOIL',           # Soil-related
        'FRLANDUSE', 'FRLAI', 'FRZ0',         # Land use categories
    }

    atmospheric_3d = []
    categorical_3d = []
    surface_2d = []
    scalar_0d = []

    for name, type_name, rank, dims, is_edge in fields:
        name_upper = name.upper()

        if rank == 3:
            # Check if this is a categorical field
            if name_upper in categorical_patterns:
                categorical_3d.append(name)
            else:
                # Assume all other 3D fields are atmospheric vertical profiles
                atmospheric_3d.append(name)
        elif rank == 2:
            # All 2D fields are surface fields
            surface_2d.append(name)
        elif rank == 0:
            # All scalar fields
            scalar_0d.append(name)

    return atmospheric_3d, categorical_3d, surface_2d, scalar_0d

def get_field_description(name):
    """
    Get a descriptive comment for a meteorological field.

    Parameters
    ----------
    name : str
        Field name

    Returns
    -------
    str
        Description string for documentation
    """
    descriptions = {
        # Common atmospheric fields
        'T': 'Temperature [K]',
        'U': 'U-wind [m/s]',
        'V': 'V-wind [m/s]',
        'QV': 'Water vapor [kg/kg]',
        'PMID': 'Mid-level pressure [Pa]',
        'AIRDEN': 'Air density [kg/m3]',
        'THETA': 'Potential temperature [K]',
        'RH': 'Relative humidity [%]',
        'OMEGA': 'Vertical velocity [Pa/s]',
        'CLDF': 'Cloud fraction [1]',
        'Z': 'Geopotential height [m]',
        'ZMID': 'Mid-level height [m]',
        'BXHEIGHT': 'Box height [m]',
        'DELP': 'Pressure thickness [Pa]',
        'AIRVOL': 'Air volume [m3]',
        'MAIRDEN': 'Moist air density [kg/m3]',
        'PEDGE': 'Edge pressure [Pa]',
        'PMID_DRY': 'Dry mid pressure [Pa]',
        'DELP_DRY': 'Dry pressure thickness [Pa]',
        'PEDGE_DRY': 'Dry edge pressure [Pa]',
        'QI': 'Ice water [kg/kg]',
        'QL': 'Liquid water [kg/kg]',
        'CMFMC': 'Cloud mass flux [kg/m2/s]',
        'DTRAIN': 'Detrainment flux [kg/m2/s]',
        'F_OF_PBL': 'Fraction in PBL [1]',
        'DQRCU': 'Conv precip rate [kg/kg/s]',
        'DQRLSAN': 'LS precip rate [kg/kg/s]',
        'PFICU': 'Conv ice precip flux [kg/m2/s]',
        'PFILSAN': 'LS ice precip flux [kg/m2/s]',
        'PFLCU': 'Conv liq precip flux [kg/m2/s]',
        'PFLLSAN': 'LS liq precip flux [kg/m2/s]',
        'TAUCLI': 'Ice cloud optical depth [1]',
        'TAUCLW': 'Water cloud optical depth [1]',
        'AIRNUMDEN': 'Air number density [molec/cm3]',
        'AVGW': 'Average molecular weight [g/mol]',
        'SPHU': 'Specific humidity [kg/kg]',
        'TV': 'Virtual temperature [K]',
        'DAIRMASS': 'Dry air mass [kg]',
        'F_UNDER_PBLTOP': 'Fraction under PBL top [1]',
        # Categorical 3D fields
        'SOILM': 'Volumetric soil moisture [m3/m3] (nsoil layers)',
        'SOILT': 'Temperature of soil layer [K] (nsoil layers)',
        'FRLANDUSE': 'Fractional land use [1] (nlanduse categories)',
        'FRSOIL': 'Fractional soil [1] (nsoil categories)',
        'FRLAI': 'LAI per land use type [m2/m2] (nlanduse categories)',
        'FRZ0': 'Roughness per land use [m] (nlanduse categories)',
        # Surface 2D fields
        'T2M': '2m temperature [K]',
        'PS': 'Surface pressure [Pa]',
        'USTAR': 'Friction velocity [m/s]',
        'U10M': '10m U-wind [m/s]',
        'V10M': '10m V-wind [m/s]',
        'PBLH': 'PBL height [m]',
        'Z0': 'Roughness length [m]',
        'LAI': 'Leaf area index [m2/m2]',
        'FRLAND': 'Land fraction [1]',
        'FROCEAN': 'Ocean fraction [1]',
        'FRSEAICE': 'Sea ice fraction [1]',
        'AREA_M2': 'Grid area [m2]',
        'SST': 'Sea surface temperature [K]',
        'TSKIN': 'Skin temperature [K]',
        'TS': 'Surface temperature [K]',
        'QV2M': '2m water vapor [kg/kg]',
        'SLP': 'Sea level pressure [Pa]',
        'PHIS': 'Surface geopotential [m2/s2]',
        'SUNCOS': 'Cosine solar zenith angle [1]',
        'SWGDN': 'Downward SW radiation [W/m2]',
        'HFLUX': 'Sensible heat flux [W/m2]',
        'EFLUX': 'Latent heat flux [W/m2]',
        'PRECCON': 'Convective precipitation [kg/m2/s]',
        'PRECLSC': 'Large-scale precipitation [kg/m2/s]',
        'PRECANV': 'Anvil precipitation [kg/m2/s]',
        'TO3': 'Total ozone [DU]',
        'TROPP': 'Tropopause pressure [Pa]',
        'TropHt': 'Tropopause height [m]'
    }

    return descriptions.get(name.upper(), f"{name} field")

def write_2d_scalar_accessor(fields, output_file):
    """
    Write a macro for 2D scalar field accessors.
    Only includes REAL type fields since the interface returns real(fp) values.

    Parameters
    ----------
    fields : list of (name, type_name, rank, dims, is_edge)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only 2D REAL fields
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 2 and type_name.upper() == 'REAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) scalar_val = this%{name}(col_i, col_j)\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")

def write_virtualmet_populate(fields, output_file):
    """
    Write a macro for populating VirtualMetType with MetState field pointers.
    Uses automatic field classification and supports all field types including scalars.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """

    # Classify fields by type and rank
    real_3d = []
    int_3d = []
    logical_3d = []
    char_3d = []
    real_2d = []
    int_2d = []
    logical_2d = []
    char_2d = []
    real_scalar = []
    int_scalar = []
    logical_scalar = []
    char_scalar = []

    for name, type_name, rank, dims, is_edge in fields:
        if rank == 3:
            if type_name.upper() == 'REAL':
                real_3d.append(name)
            elif type_name.upper() == 'INTEGER':
                int_3d.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_3d.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_3d.append(name)
        elif rank == 2:
            if type_name.upper() == 'REAL':
                real_2d.append(name)
            elif type_name.upper() == 'INTEGER':
                int_2d.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_2d.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_2d.append(name)
        elif rank == 0:  # Include scalar fields
            if type_name.upper() == 'REAL':
                real_scalar.append(name)
            elif type_name.upper() == 'INTEGER':
                int_scalar.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_scalar.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_scalar.append(name)

    with open(output_file, 'w') as f:
        f.write("! Generated macro for populating VirtualMetType with MetState field pointers\n")
        f.write("! This macro should be included in the populate_virtual_column subroutine\n")
        f.write("! Auto-generated from MetState field definitions with native type support\n\n")

        # Populate 3D REAL field pointers
        if real_3d:
            f.write("! Populate 3D REAL atmospheric field pointers (vertical levels: nlev)\n")
            for name in sorted(real_3d):
                f.write(f"call this%met_state%get_field_ptr('{name}', grid_i, grid_j, &\n")
                f.write(f"                                 col_ptr=column_ptr, rc=field_rc)\n")
                f.write(f"if (field_rc == 0 .and. associated(column_ptr)) then\n")
                f.write(f"   virtual_col%met%{name} => column_ptr\n")
                f.write(f"end if\n\n")

        # Populate 3D INTEGER field pointers
        if int_3d:
            f.write("! Populate 3D INTEGER field pointers\n")
            for name in sorted(int_3d):
                f.write(f"call this%met_state%get_field_ptr_int('{name}', grid_i, grid_j, &\n")
                f.write(f"                                      col_ptr=column_ptr_int, rc=field_rc)\n")
                f.write(f"if (field_rc == 0 .and. associated(column_ptr_int)) then\n")
                f.write(f"   virtual_col%met%{name} => column_ptr_int\n")
                f.write(f"end if\n\n")

        # Populate 3D LOGICAL field pointers
        if logical_3d:
            f.write("! Populate 3D LOGICAL field pointers\n")
            for name in sorted(logical_3d):
                f.write(f"call this%met_state%get_field_ptr_logical('{name}', grid_i, grid_j, &\n")
                f.write(f"                                          col_ptr=column_ptr_logical, rc=field_rc)\n")
                f.write(f"if (field_rc == 0 .and. associated(column_ptr_logical)) then\n")
                f.write(f"   virtual_col%met%{name} => column_ptr_logical\n")
                f.write(f"end if\n\n")

        # Populate 3D CHARACTER field pointers
        if char_3d:
            f.write("! Populate 3D CHARACTER field pointers\n")
            for name in sorted(char_3d):
                f.write(f"! Note: CHARACTER field {name} accessed directly from MetState\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   virtual_col%met%{name} => this%met_state%{name}(:)\n")
                f.write(f"end if\n\n")

        # Populate 2D REAL scalar fields
        if real_2d:
            f.write("! Populate 2D REAL scalar fields\n")
            for name in sorted(real_2d):
                f.write(f"call this%met_state%get_field_ptr('{name}', grid_i, grid_j, &\n")
                f.write(f"                                 scalar_val=scalar_val, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val\n")
                f.write(f"end if\n\n")

        # Populate 2D INTEGER scalar fields
        if int_2d:
            f.write("! Populate 2D INTEGER scalar fields\n")
            for name in sorted(int_2d):
                f.write(f"call this%met_state%get_field_ptr_int('{name}', grid_i, grid_j, &\n")
                f.write(f"                                      scalar_val=scalar_val_int, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val_int\n")
                f.write(f"end if\n\n")

        # Populate 2D LOGICAL scalar fields
        if logical_2d:
            f.write("! Populate 2D LOGICAL scalar fields\n")
            for name in sorted(logical_2d):
                f.write(f"call this%met_state%get_field_ptr_logical('{name}', grid_i, grid_j, &\n")
                f.write(f"                                          scalar_val=scalar_val_logical, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val_logical\n")
                f.write(f"end if\n\n")

        # Populate 2D CHARACTER scalar fields
        if char_2d:
            f.write("! Populate 2D CHARACTER scalar fields\n")
            for name in sorted(char_2d):
                f.write(f"! Note: CHARACTER field {name} accessed directly from MetState\n")
                f.write(f"virtual_col%met%{name} = this%met_state%{name}(grid_i, grid_j)\n\n")

        # Populate scalar REAL fields
        if real_scalar:
            f.write("! Populate scalar REAL fields\n")
            for name in sorted(real_scalar):
                f.write(f"call this%met_state%get_field_ptr('{name}', &\n")
                f.write(f"                                 scalar_val=scalar_val, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val\n")
                f.write(f"end if\n\n")

        # Populate scalar INTEGER fields
        if int_scalar:
            f.write("! Populate scalar INTEGER fields\n")
            for name in sorted(int_scalar):
                f.write(f"call this%met_state%get_field_ptr_int('{name}', &\n")
                f.write(f"                                      scalar_val=scalar_val_int, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val_int\n")
                f.write(f"end if\n\n")

        # Populate scalar LOGICAL fields
        if logical_scalar:
            f.write("! Populate scalar LOGICAL fields\n")
            for name in sorted(logical_scalar):
                f.write(f"call this%met_state%get_field_ptr_logical('{name}', &\n")
                f.write(f"                                          scalar_val=scalar_val_logical, rc=field_rc)\n")
                f.write(f"if (field_rc == 0) then\n")
                f.write(f"   virtual_col%met%{name} = scalar_val_logical\n")
                f.write(f"end if\n\n")

        # Populate scalar CHARACTER fields
        if char_scalar:
            f.write("! Populate scalar CHARACTER fields\n")
            for name in sorted(char_scalar):
                f.write(f"! Note: CHARACTER field {name} accessed directly from MetState\n")
                f.write(f"virtual_col%met%{name} = this%met_state%{name}\n\n")

        f.write("! Note: All field types (REAL, INTEGER, LOGICAL, CHARACTER) including scalars are now supported.\n")
        f.write("! Type-specific accessor functions used for numeric/logical fields, direct access for CHARACTER fields.\n")

def write_virtualmet_type(fields, output_file):
    """
    Write the VirtualMetType definition macro based on MetState field definitions.
    Uses automatic field classification and supports all field types including scalars.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """

    # Automatically classify fields by type and rank
    atmospheric_3d, categorical_3d, surface_2d, scalar_0d = classify_fields(fields)

    # Further classify by type and rank
    real_3d = []
    int_3d = []
    logical_3d = []
    char_3d = []
    real_2d = []
    int_2d = []
    logical_2d = []
    char_2d = []
    real_scalar = []
    int_scalar = []
    logical_scalar = []
    char_scalar = []

    for name, type_name, rank, dims, is_edge in fields:
        if rank == 3:
            if type_name.upper() == 'REAL':
                real_3d.append(name)
            elif type_name.upper() == 'INTEGER':
                int_3d.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_3d.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_3d.append(name)
        elif rank == 2:
            if type_name.upper() == 'REAL':
                real_2d.append(name)
            elif type_name.upper() == 'INTEGER':
                int_2d.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_2d.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_2d.append(name)
        elif rank == 0:  # Include scalar fields now
            if type_name.upper() == 'REAL':
                real_scalar.append(name)
            elif type_name.upper() == 'INTEGER':
                int_scalar.append(name)
            elif type_name.upper() == 'LOGICAL':
                logical_scalar.append(name)
            elif type_name.upper() == 'CHARACTER':
                char_scalar.append(name)

    with open(output_file, 'w') as f:
        f.write("! Generated VirtualMetType definition based on MetState field definitions\n")
        f.write("! This macro should be included in the VirtualColumn_Mod.F90 type definition\n")
        f.write("! Auto-generated from MetState field definitions with full field type support\n\n")

        f.write("   !> \\brief Virtual meteorological data container with direct pointers and scalar values\n")
        f.write("   !! \\details Contains pointers to meteorological fields for a single column plus scalar values.\n")
        f.write("   !! Pointers are set directly from MetState accessor functions, eliminating\n")
        f.write("   !! data copying and providing efficient field access.\n")
        f.write("   type :: VirtualMetType\n")

        # Write 3D REAL field pointers (vertical profiles)
        if real_3d:
            f.write("      ! 3D REAL atmospheric fields (vertical profiles: nlev) - pointers to MetState data\n")
            for name in sorted(real_3d):
                comment = get_field_description(name)
                f.write(f"      real(fp), pointer :: {name}(:) => null()  !< {comment}\n")
            f.write("\n")

        # Write 3D INTEGER field pointers
        if int_3d:
            f.write("      ! 3D INTEGER fields - pointers to MetState data\n")
            for name in sorted(int_3d):
                comment = get_field_description(name)
                f.write(f"      integer, pointer :: {name}(:) => null()  !< {comment}\n")
            f.write("\n")

        # Write 3D LOGICAL field pointers
        if logical_3d:
            f.write("      ! 3D LOGICAL fields - pointers to MetState data\n")
            for name in sorted(logical_3d):
                comment = get_field_description(name)
                f.write(f"      logical, pointer :: {name}(:) => null()  !< {comment}\n")
            f.write("\n")

        # Write 3D CHARACTER field pointers
        if char_3d:
            f.write("      ! 3D CHARACTER fields - pointers to MetState data\n")
            for name in sorted(char_3d):
                comment = get_field_description(name)
                f.write(f"      character(len=255), pointer :: {name}(:) => null()  !< {comment}\n")
            f.write("\n")

        # Write 2D REAL surface fields (scalars)
        if real_2d:
            f.write("      ! 2D REAL surface fields (scalars) - direct values from MetState\n")
            for name in sorted(real_2d):
                comment = get_field_description(name)
                f.write(f"      real(fp) :: {name}  !< {comment}\n")
            f.write("\n")

        # Write 2D INTEGER surface fields
        if int_2d:
            f.write("      ! 2D INTEGER surface fields - direct values from MetState\n")
            for name in sorted(int_2d):
                comment = get_field_description(name)
                f.write(f"      integer :: {name}  !< {comment}\n")
            f.write("\n")

        # Write 2D LOGICAL surface fields
        if logical_2d:
            f.write("      ! 2D LOGICAL surface fields - direct values from MetState\n")
            for name in sorted(logical_2d):
                comment = get_field_description(name)
                f.write(f"      logical :: {name}  !< {comment}\n")
            f.write("\n")

        # Write 2D CHARACTER surface fields
        if char_2d:
            f.write("      ! 2D CHARACTER surface fields - direct values from MetState\n")
            for name in sorted(char_2d):
                comment = get_field_description(name)
                f.write(f"      character(len=255) :: {name}  !< {comment}\n")
            f.write("\n")

        # Write scalar REAL fields
        if real_scalar:
            f.write("      ! Scalar REAL fields - direct values from MetState\n")
            for name in sorted(real_scalar):
                comment = get_field_description(name)
                f.write(f"      real(fp) :: {name}  !< {comment}\n")
            f.write("\n")

        # Write scalar INTEGER fields
        if int_scalar:
            f.write("      ! Scalar INTEGER fields - direct values from MetState\n")
            for name in sorted(int_scalar):
                comment = get_field_description(name)
                f.write(f"      integer :: {name}  !< {comment}\n")
            f.write("\n")

        # Write scalar LOGICAL fields
        if logical_scalar:
            f.write("      ! Scalar LOGICAL fields - direct values from MetState\n")
            for name in sorted(logical_scalar):
                comment = get_field_description(name)
                f.write(f"      logical :: {name}  !< {comment}\n")
            f.write("\n")

        # Write scalar CHARACTER fields
        if char_scalar:
            f.write("      ! Scalar CHARACTER fields - direct values from MetState\n")
            for name in sorted(char_scalar):
                comment = get_field_description(name)
                f.write(f"      character(len=255) :: {name}  !< {comment}\n")
            f.write("\n")

        f.write("   contains\n")
        f.write("      procedure :: cleanup => virtual_met_cleanup\n")
        f.write("   end type VirtualMetType\n\n")

def write_virtualmet_cleanup(fields, output_file):
    """
    Write the VirtualMetType cleanup procedure macro.
    Uses automatic field classification instead of hardcoded lists.
    Now includes scalar fields since VirtualMet contains them.

    Parameters
    ----------
    fields : list of tuple
        List of (name, type_name, rank, dims, is_edge) for each field.
    output_file : str
        Path to output .inc file.
    """

    # Automatically classify fields
    atmospheric_3d, categorical_3d, surface_2d, scalar_0d = classify_fields(fields)

    # Group scalar fields by type for proper initialization
    scalar_real = []
    scalar_int = []
    scalar_logical = []

    for name, type_name, rank, dims, is_edge in fields:
        if name in scalar_0d:
            if type_name == 'real':
                scalar_real.append(name)
            elif type_name == 'integer':
                scalar_int.append(name)
            elif type_name == 'logical':
                scalar_logical.append(name)

    with open(output_file, 'w') as f:
        f.write("! Generated VirtualMetType cleanup procedure\n")
        f.write("! This macro should be included in the virtual_met_cleanup subroutine\n")
        f.write("! Auto-generated from MetState field definitions\n")
        f.write("! Now includes scalar fields since VirtualMet contains them\n")
        f.write("! Now includes scalar fields since VirtualMet contains them\n\n")

        # Nullify 3D atmospheric field pointers
        if atmospheric_3d:
            f.write("      ! Nullify 3D atmospheric field pointers (do not deallocate - they point to MetState data)\n")
            for name in sorted(atmospheric_3d):
                f.write(f"      this%{name} => null()\n")

        # Nullify categorical 3D field pointers
        if categorical_3d:
            f.write("\n      ! Nullify 3D categorical field pointers (do not deallocate - they point to MetState data)\n")
            for name in sorted(categorical_3d):
                f.write(f"      this%{name} => null()\n")

        # Reset surface 2D fields by type
        if surface_2d:
            # Separate surface fields by type
            surface_real = []
            surface_int = []
            surface_logical = []

            for name, type_name, rank, dims, is_edge in fields:
                if name in surface_2d:
                    if type_name == 'real':
                        surface_real.append(name)
                    elif type_name == 'integer':
                        surface_int.append(name)
                    elif type_name == 'logical':
                        surface_logical.append(name)

            if surface_real:
                f.write("\n      ! Reset 2D real scalar fields to default values\n")
                for name in sorted(surface_real):
                    f.write(f"      this%{name} = 0.0_fp\n")

            if surface_int:
                f.write("\n      ! Reset 2D integer scalar fields to default values\n")
                for name in sorted(surface_int):
                    f.write(f"      this%{name} = 0\n")

            if surface_logical:
                f.write("\n      ! Reset 2D logical scalar fields to default values\n")
                for name in sorted(surface_logical):
                    f.write(f"      this%{name} = .false.\n")

        # Reset scalar fields by type (now included since VirtualMet contains them)
        if scalar_real:
            f.write("\n      ! Reset real scalar fields to default values\n")
            for name in sorted(scalar_real):
                f.write(f"      this%{name} = 0.0_fp\n")

        if scalar_int:
            f.write("\n      ! Reset integer scalar fields to default values\n")
            for name in sorted(scalar_int):
                f.write(f"      this%{name} = 0\n")

        if scalar_logical:
            f.write("\n      ! Reset logical scalar fields to default values\n")
            for name in sorted(scalar_logical):
                f.write(f"      this%{name} = .false.\n")

        f.write("\n")

def write_virtualcolumn_populate(fields, output_file):
    """
    Write a macro for populating VirtualColumn with all MetState fields (old approach).
    Uses automatic field classification instead of hardcoded lists.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """

    # Automatically classify fields
    atmospheric_3d, categorical_3d, surface_2d, scalar_0d = classify_fields(fields)

    with open(output_file, 'w') as f:
        f.write("! Generated macro for populating VirtualColumn with MetState fields\n")
        f.write("! This macro should be included in the populate_virtual_column subroutine\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Add 3D vertical profile fields only (exclude categorical 3D fields)
        if atmospheric_3d:
            f.write("! Add 3D atmospheric vertical profile fields\n")
            for name in sorted(atmospheric_3d):
                f.write(f"! Add {name} 3D vertical profile field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   call virtual_col%add_met_field('{name}', this%met_state%{name}(grid_i, grid_j, :), rc)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")

        # Add 2D surface fields
        if surface_2d:
            f.write("! Add 2D surface fields\n")
            for name in sorted(surface_2d):
                f.write(f"! Add {name} 2D surface field\n")
                f.write(f"if (allocated(this%met_state%{name})) then\n")
                f.write(f"   allocate(temp_column(1))\n")
                f.write(f"   temp_column(1) = this%met_state%{name}(grid_i, grid_j)\n")
                f.write(f"   call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"   deallocate(temp_column)\n")
                f.write(f"   if (rc /= CC_SUCCESS) return\n")
                f.write(f"end if\n\n")

        # Add scalar fields
        if scalar_0d:
            f.write("! Add scalar fields\n")
            for name in sorted(scalar_0d):
                f.write(f"! Add {name} scalar field\n")
                f.write(f"allocate(temp_column(1))\n")
                f.write(f"temp_column(1) = this%met_state%{name}\n")
                f.write(f"call virtual_col%add_met_field('{name}', temp_column, rc)\n")
                f.write(f"deallocate(temp_column)\n")
                f.write(f"if (rc /= CC_SUCCESS) return\n\n")

        f.write("! Note: Categorical 3D fields are automatically detected and excluded.\n")
        f.write("! They require specific category indices, not column extraction.\n")
        f.write("! Processes needing these should access them directly from MetState.\n")

def write_set_field_2d_real(fields, output_file):
    """
    Write a macro for setting 2D REAL field values.
    Now handles case where arrays may not be allocated.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 2D REAL MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n")
        f.write("! Now handles arrays that may not be allocated\n\n")

        # Generate cases for 2D REAL fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 2 and type_name == 'real':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
                # Check if field needs special allocation parameters
                if name.lower() in ("frlanduse", "frlai", "frz0"):
                    f.write(f"   if (.not. allocated(this%{name})) then\n")
                    f.write(f"      ! Allocate with default NSURFTYPE if not already allocated\n")
                    f.write(f"      if (this%NSURFTYPE > 0) then\n")
                    f.write(f"         call this%allocate_arrays('{name}', error_mgr, rc)\n")
                    f.write(f"         if (rc /= CC_SUCCESS) return\n")
                    f.write(f"      else\n")
                    f.write(f"         call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                    f.write(f"            'Cannot allocate {name}: NSURFTYPE not set', rc)\n")
                    f.write(f"         return\n")
                    f.write(f"      end if\n")
                    f.write(f"   end if\n")
                elif name.lower() in ("frsoil", "soilm", "soilt"):
                    f.write(f"   if (.not. allocated(this%{name})) then\n")
                    f.write(f"      ! Allocate with default soil parameters if not already allocated\n")
                    f.write(f"      if (this%nSOILTYPE > 0) then\n")
                    f.write(f"         call this%allocate_arrays('{name}', error_mgr, rc)\n")
                    f.write(f"         if (rc /= CC_SUCCESS) return\n")
                    f.write(f"      else\n")
                    f.write(f"         call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                    f.write(f"            'Cannot allocate {name}: soil parameters not set', rc)\n")
                    f.write(f"         return\n")
                    f.write(f"      end if\n")
                    f.write(f"   end if\n")
                else:
                    f.write(f"   if (.not. allocated(this%{name})) then\n")
                    f.write(f"      call this%allocate_arrays('{name}', error_mgr, rc)\n")
                    f.write(f"      if (rc /= CC_SUCCESS) return\n")
                    f.write(f"   end if\n")
                f.write(f"   this%{name} = field_data\n")
                f.write(f"   rc = CC_SUCCESS\n\n")

def write_set_field_2d_int(fields, output_file):
    """
    Write a macro for setting 2D INTEGER field values.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 2D INTEGER MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Generate cases for 2D INTEGER fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 2 and type_name == 'integer':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
                f.write(f"   if (.not. allocated(this%{name})) then\n")
                f.write(f"      call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                f.write(f"         'Field {name} not allocated', rc)\n")
                f.write(f"      return\n")
                f.write(f"   end if\n")
                f.write(f"   this%{name} = field_data\n")
                f.write(f"   rc = CC_SUCCESS\n\n")

def write_set_field_2d_logical(fields, output_file):
    """
    Write a macro for setting 2D LOGICAL field values.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 2D LOGICAL MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Generate cases for 2D LOGICAL fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 2 and type_name == 'logical':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
                f.write(f"   if (.not. allocated(this%{name})) then\n")
                f.write(f"      call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                f.write(f"         'Field {name} not allocated', rc)\n")
                f.write(f"      return\n")
                f.write(f"   end if\n")
                f.write(f"   this%{name} = field_data\n")
                f.write(f"   rc = CC_SUCCESS\n\n")

def get_conditional_allocation_info(field_name):
    """
    Determine if a field requires conditional allocation and return allocation details.

    Parameters
    ----------
    field_name : str
        Name of the field to check

    Returns
    -------
    dict or None
        Dictionary with 'condition', 'dimension', and 'type' if conditional,
        None if standard allocation
    """
    field_lower = field_name.lower()

    # Define conditional allocation patterns
    conditional_patterns = {
        'soilm': {
            'condition': 'nsoil > 0',
            'dimension': 'nsoil',
            'dimension_var': 'this%nSOIL',
            'type': 'soil'
        },
        'soilt': {
            'condition': 'nsoil > 0',
            'dimension': 'nsoil',
            'dimension_var': 'this%nSOIL',
            'type': 'soil'
        },
        'frsoil': {
            'condition': 'nsoiltype > 0',
            'dimension': 'nsoiltype',
            'dimension_var': 'this%nSOILTYPE',
            'type': 'soiltype'
        },
        'iland': {
            'condition': 'nSURFTYPE > 0',
            'dimension': 'nSURFTYPE',
            'dimension_var': 'this%NSURFTYPE',
            'type': 'surftype'
        }
    }

    # Check for surface type fields (fields starting with 'fr' that aren't soil-related)
    if field_lower.startswith('fr') and field_lower not in ['frsoil']:
        # Common surface type fields: frlanduse, frlai, frz0, etc.
        conditional_patterns[field_lower] = {
            'condition': 'nSURFTYPE > 0',
            'dimension': 'nSURFTYPE',
            'dimension_var': 'this%NSURFTYPE',
            'type': 'surftype'
        }

    return conditional_patterns.get(field_lower)

def write_set_field_3d_real(fields, output_file):
    """
    Write a macro for setting 3D REAL field values.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 3D REAL MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Generate cases for 3D REAL fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 3 and type_name == 'real':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")

                # Check if this field has conditional allocation
                conditional_info = get_conditional_allocation_info(name)

                if conditional_info:
                    # Use automatic allocation on assignment for conditional fields
                    f.write(f"   ! Automatic allocation on assignment for conditional field {name}\n")
                    f.write(f"   ! Field type: {conditional_info['type']}\n")
                    f.write(f"   this%{name} = field_data\n")
                    f.write(f"   rc = CC_SUCCESS\n\n")
                else:
                    # Standard allocation check for other fields
                    f.write(f"   if (.not. allocated(this%{name})) then\n")
                    f.write(f"      call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                    f.write(f"         'Field {name} not allocated', rc)\n")
                    f.write(f"      return\n")
                    f.write(f"   end if\n")
                    f.write(f"   this%{name} = field_data\n")
                    f.write(f"   rc = CC_SUCCESS\n\n")

def write_set_field_3d_int(fields, output_file):
    """
    Write a macro for setting 3D INTEGER field values.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 3D INTEGER MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Generate cases for 3D INTEGER fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 3 and type_name == 'integer':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
                f.write(f"   if (.not. allocated(this%{name})) then\n")
                f.write(f"      call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                f.write(f"         'Field {name} not allocated', rc)\n")
                f.write(f"      return\n")
                f.write(f"   end if\n")
                f.write(f"   this%{name} = field_data\n")
                f.write(f"   rc = CC_SUCCESS\n\n")

def write_set_field_3d_logical(fields, output_file):
    """
    Write a macro for setting 3D LOGICAL field values.
    """
    with open(output_file, 'w') as f:
        f.write("! Generated macro for setting 3D LOGICAL MetState field values\n")
        f.write("! Auto-generated from MetState field definitions\n\n")

        # Generate cases for 3D LOGICAL fields
        for name, type_name, rank, dims, is_edge in fields:
            if rank == 3 and type_name == 'logical':
                labels = sorted({name, name.lower()})
                f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
                f.write(f"   if (.not. allocated(this%{name})) then\n")
                f.write(f"      call error_mgr%report_error(ERROR_INVALID_INPUT, &\n")
                f.write(f"         'Field {name} not allocated', rc)\n")
                f.write(f"      return\n")
                f.write(f"   end if\n")
                f.write(f"   this%{name} = field_data\n")
                f.write(f"   rc = CC_SUCCESS\n\n")

def write_scalar_accessor(fields, output_file):
    """
    Write a macro for scalar (rank-0) field accessors.
    Only includes REAL type fields since the interface returns real(fp) values.

    Parameters
    ----------
    fields : list of (name, type_name, rank, dims, is_edge)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'REAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
            f.write(f"   rc = 0\n")
            f.write(f"   return\n")

def write_column_accessor_int(fields, output_file):
    """
    Write a macro for 3D column field accessors for INTEGER fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 3 and type_name.upper() == 'INTEGER':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) column_ptr => this%{name}(col_i, col_j, :)\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_column_accessor_logical(fields, output_file):
    """
    Write a macro for 3D column field accessors for LOGICAL fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 3 and type_name.upper() == 'LOGICAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) column_ptr => this%{name}(col_i, col_j, :)\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_2d_scalar_accessor_int(fields, output_file):
    """
    Write a macro for 2D scalar field accessors for INTEGER fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 2 and type_name.upper() == 'INTEGER':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) scalar_val = this%{name}(col_i, col_j)\n")
        f.write("case default\n")
        f.write("   scalar_val = 0\n")

def write_2d_scalar_accessor_logical(fields, output_file):
    """
    Write a macro for 2D scalar field accessors for LOGICAL fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 2 and type_name.upper() == 'LOGICAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) scalar_val = this%{name}(col_i, col_j)\n")
        f.write("case default\n")
        f.write("   scalar_val = .false.\n")

def write_scalar_accessor_int(fields, output_file):
    """
    Write a macro for scalar (rank-0) field accessors for INTEGER fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'INTEGER':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
        f.write("case default\n")
        f.write("   scalar_val = 0\n")

def write_scalar_accessor_logical(fields, output_file):
    """
    Write a macro for scalar (rank-0) field accessors for LOGICAL fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'LOGICAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
        f.write("case default\n")
        f.write("   scalar_val = .false.\n")

def write_set_field_scalar_real(fields, output_file):
    """
    Write a macro for setting scalar REAL fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'REAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   this%{name} = field_data\n")

def write_set_field_scalar_int(fields, output_file):
    """
    Write a macro for setting scalar INTEGER fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'INTEGER':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   this%{name} = field_data\n")

def write_set_field_scalar_logical(fields, output_file):
    """
    Write a macro for setting scalar LOGICAL fields.
    """
    unique_fields = {}
    for name, type_name, rank, dims, is_edge in fields:
        if rank == 0 and type_name.upper() == 'LOGICAL':
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, type_name, rank, dims, is_edge)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, type_name, rank, dims, is_edge = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   this%{name} = field_data\n")

def write_multiple_fields_interface(fields, output_file):
    """
    Generate a comprehensive set_multiple_fields function with all fields as optional arguments.
    This allows setting multiple MetState fields in a single call directly from host model data.
    """

    # Group fields by type and rank
    real_2d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                      if type_name == 'real' and rank == 2]
    real_3d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                      if type_name == 'real' and rank == 3]
    int_2d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                     if type_name == 'integer' and rank == 2]
    int_3d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                     if type_name == 'integer' and rank == 3]
    logical_2d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                         if type_name == 'logical' and rank == 2]
    logical_3d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                         if type_name == 'logical' and rank == 3]
    char_2d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                      if type_name == 'character' and rank == 2]
    char_3d_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                      if type_name == 'character' and rank == 3]
    real_scalar_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                          if type_name == 'real' and rank == 0]
    int_scalar_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                         if type_name == 'integer' and rank == 0]
    logical_scalar_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                             if type_name == 'logical' and rank == 0]
    char_scalar_fields = [(name, rank, dims) for name, type_name, rank, dims, is_edge in fields
                          if type_name == 'character' and rank == 0]

    with open(output_file, 'w') as f:
        f.write("!> \\brief Set multiple meteorological fields directly from host model data\n")
        f.write("!!\n")
        f.write("!! This subroutine allows setting multiple MetState fields in a single call\n")
        f.write("!! with data directly from host model arrays. Only provide the optional\n")
        f.write("!! arguments for fields you want to set. This is the most efficient way\n")
        f.write("!! to transfer data from host models to CATChem MetState.\n")
        f.write("!!\n")
        f.write("!! \\param[inout] this      MetStateType object\n")
        f.write("!! \\param[in]    field_names Array of field names to set\n")
        f.write("!! \\param[inout] error_mgr Error manager for reporting\n")
        f.write("!! \\param[out]   rc        Return code\n")

        # Write subroutine header
        f.write("subroutine metstate_set_multiple_fields(this, field_names, error_mgr, rc")

        # Count total optional arguments
        all_field_groups = [real_2d_fields, int_2d_fields, logical_2d_fields, char_2d_fields,
                           real_3d_fields, int_3d_fields, logical_3d_fields, char_3d_fields,
                           real_scalar_fields, int_scalar_fields, logical_scalar_fields, char_scalar_fields]
        has_fields = any(len(group) > 0 for group in all_field_groups)

        if has_fields:
            f.write(", &\n")

            # Write all the optional field arguments with proper continuation
            groups_written = []

            # 2D REAL fields
            if real_2d_fields:
                f.write("                                        ! 2D REAL fields\n")
                for i, (field_name, rank, dims) in enumerate(real_2d_fields):
                    if i == len(real_2d_fields) - 1 and not (int_2d_fields or logical_2d_fields or real_3d_fields or int_3d_fields or logical_3d_fields or real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('real_2d')

            # 2D INTEGER fields
            if int_2d_fields:
                f.write("                                        ! 2D INTEGER fields\n")
                for i, (field_name, rank, dims) in enumerate(int_2d_fields):
                    if i == len(int_2d_fields) - 1 and not (logical_2d_fields or real_3d_fields or int_3d_fields or logical_3d_fields or real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('int_2d')

            # 2D LOGICAL fields
            if logical_2d_fields:
                f.write("                                        ! 2D LOGICAL fields\n")
                for i, (field_name, rank, dims) in enumerate(logical_2d_fields):
                    if i == len(logical_2d_fields) - 1 and not (real_3d_fields or int_3d_fields or logical_3d_fields or real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('logical_2d')

            # 3D REAL fields
            if real_3d_fields:
                f.write("                                        ! 3D REAL fields\n")
                for i, (field_name, rank, dims) in enumerate(real_3d_fields):
                    if i == len(real_3d_fields) - 1 and not (int_3d_fields or logical_3d_fields or real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('real_3d')

            # 3D INTEGER fields
            if int_3d_fields:
                f.write("                                        ! 3D INTEGER fields\n")
                for i, (field_name, rank, dims) in enumerate(int_3d_fields):
                    if i == len(int_3d_fields) - 1 and not (logical_3d_fields or real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('int_3d')

            # 3D LOGICAL fields
            if logical_3d_fields:
                f.write("                                        ! 3D LOGICAL fields\n")
                for i, (field_name, rank, dims) in enumerate(logical_3d_fields):
                    if i == len(logical_3d_fields) - 1 and not (real_scalar_fields or int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('logical_3d')

            # Scalar REAL fields
            if real_scalar_fields:
                f.write("                                        ! Scalar REAL fields\n")
                for i, (field_name, rank, dims) in enumerate(real_scalar_fields):
                    if i == len(real_scalar_fields) - 1 and not (int_scalar_fields or logical_scalar_fields):
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('real_scalar')

            # Scalar INTEGER fields
            if int_scalar_fields:
                f.write("                                        ! Scalar INTEGER fields\n")
                for i, (field_name, rank, dims) in enumerate(int_scalar_fields):
                    if i == len(int_scalar_fields) - 1 and not logical_scalar_fields:
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('int_scalar')

            # Scalar LOGICAL fields
            if logical_scalar_fields:
                f.write("                                        ! Scalar LOGICAL fields\n")
                for i, (field_name, rank, dims) in enumerate(logical_scalar_fields):
                    if i == len(logical_scalar_fields) - 1:
                        f.write(f"                                        {field_name}_data)\n")
                    else:
                        f.write(f"                                        {field_name}_data, &\n")
                groups_written.append('logical_scalar')
        else:
            f.write(")\n")

        # Write the argument declarations
        f.write("   use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_NOT_FOUND\n")
        f.write("   implicit none\n")
        f.write("   class(MetStateType), intent(inout) :: this\n")
        f.write("   character(len=*), intent(in) :: field_names(:)\n")
        f.write("   type(ErrorManagerType), pointer, intent(inout) :: error_mgr\n")
        f.write("   integer, intent(out) :: rc\n\n")

        # Write all the optional argument declarations
        if real_2d_fields:
            f.write("   ! Optional 2D REAL field arguments\n")
            for field_name, rank, dims in real_2d_fields:
                f.write(f"   real(fp), intent(in), optional :: {field_name}_data(:,:)\n")
            f.write("\n")

        if int_2d_fields:
            f.write("   ! Optional 2D INTEGER field arguments\n")
            for field_name, rank, dims in int_2d_fields:
                f.write(f"   integer, intent(in), optional :: {field_name}_data(:,:)\n")
            f.write("\n")

        if logical_2d_fields:
            f.write("   ! Optional 2D LOGICAL field arguments\n")
            for field_name, rank, dims in logical_2d_fields:
                f.write(f"   logical, intent(in), optional :: {field_name}_data(:,:)\n")
            f.write("\n")

        if real_3d_fields:
            f.write("   ! Optional 3D REAL field arguments\n")
            for field_name, rank, dims in real_3d_fields:
                f.write(f"   real(fp), intent(in), optional :: {field_name}_data(:,:,:)\n")
            f.write("\n")

        if int_3d_fields:
            f.write("   ! Optional 3D INTEGER field arguments\n")
            for field_name, rank, dims in int_3d_fields:
                f.write(f"   integer, intent(in), optional :: {field_name}_data(:,:,:)\n")
            f.write("\n")

        if logical_3d_fields:
            f.write("   ! Optional 3D LOGICAL field arguments\n")
            for field_name, rank, dims in logical_3d_fields:
                f.write(f"   logical, intent(in), optional :: {field_name}_data(:,:,:)\n")
            f.write("\n")

        if real_scalar_fields:
            f.write("   ! Optional scalar REAL field arguments\n")
            for field_name, rank, dims in real_scalar_fields:
                f.write(f"   real(fp), intent(in), optional :: {field_name}_data\n")
            f.write("\n")

        if int_scalar_fields:
            f.write("   ! Optional scalar INTEGER field arguments\n")
            for field_name, rank, dims in int_scalar_fields:
                f.write(f"   integer, intent(in), optional :: {field_name}_data\n")
            f.write("\n")

        if logical_scalar_fields:
            f.write("   ! Optional scalar LOGICAL field arguments\n")
            for field_name, rank, dims in logical_scalar_fields:
                f.write(f"   logical, intent(in), optional :: {field_name}_data\n")
            f.write("\n")

        # Write local variables
        f.write("   ! Local variables\n")
        f.write("   integer :: i, local_rc\n")
        f.write("   character(len=64) :: current_field\n\n")

        # Write the main logic
        f.write("   rc = CC_SUCCESS\n\n")
        f.write("   ! Loop through requested field names and set corresponding data\n")
        f.write("   do i = 1, size(field_names)\n")
        f.write("      current_field = trim(adjustl(field_names(i)))\n")
        f.write("      local_rc = CC_SUCCESS\n\n")
        f.write("      select case (current_field)\n")

        # Generate case statements for all fields
        all_fields = (real_2d_fields + int_2d_fields + logical_2d_fields +
                     real_3d_fields + int_3d_fields + logical_3d_fields +
                     real_scalar_fields + int_scalar_fields + logical_scalar_fields)

        for field_name, rank, dims in all_fields:
            f.write(f"      case ('{field_name}', '{field_name.lower()}')\n")
            f.write(f"         if (present({field_name}_data)) then\n")
            f.write(f"            call this%set_field('{field_name}', {field_name}_data, error_mgr, local_rc)\n")
            f.write(f"         else\n")
            f.write(f"            call error_mgr%report_error(ERROR_NOT_FOUND, &\n")
            f.write(f"               'Field {field_name} requested but data not provided', local_rc)\n")
            f.write(f"            local_rc = CC_FAILURE\n")
            f.write(f"         endif\n")

        f.write("      case default\n")
        f.write("         call error_mgr%report_error(ERROR_NOT_FOUND, &\n")
        f.write("            'Unknown field name: ' // trim(current_field), local_rc)\n")
        f.write("         local_rc = CC_FAILURE\n")
        f.write("      end select\n\n")

        f.write("      if (local_rc /= CC_SUCCESS) then\n")
        f.write("         write(*,'(A,A)') 'Warning: Failed to set field: ', trim(current_field)\n")
        f.write("         rc = CC_FAILURE\n")
        f.write("         ! Continue with other fields\n")
        f.write("      endif\n")
        f.write("   end do\n\n")

        f.write("end subroutine metstate_set_multiple_fields\n")

def main():
    """
    Main entry point for the macro generator script.

    Parses command-line arguments and dispatches to the appropriate macro writer.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate Fortran macros for MetStateType fields.",
        epilog="Example: generate_metstate_macros.py allocate metstate_mod.F90 metstate_allocate_fields.inc"
    )
    parser.add_argument('mode', choices=['accessor', 'allocate', 'deallocate', 'column_accessor', '2d_scalar_accessor', 'scalar_accessor', 'column_accessor_int', '2d_scalar_accessor_int', 'scalar_accessor_int', 'column_accessor_logical', '2d_scalar_accessor_logical', 'scalar_accessor_logical', 'set_field_2d_real', 'set_field_2d_int', 'set_field_2d_logical', 'set_field_3d_real', 'set_field_3d_int', 'set_field_3d_logical', 'set_field_scalar_real', 'set_field_scalar_int', 'set_field_scalar_logical', 'multiple_fields_interface', 'virtualcolumn_populate', 'virtualmet_populate', 'virtualmet_type', 'virtualmet_cleanup'],
                        help="Type of macro to generate.")
    parser.add_argument('input_f90', help="Input Fortran file containing MetStateType.")
    parser.add_argument('output_inc', help="Output .inc file to write.")
    args = parser.parse_args()

    fields = parse_metstate_type(args.input_f90)
    if not fields:
        print(f"Warning: No fields found in {args.input_f90}.")

    if args.mode == 'accessor':
        write_accessor(fields, args.output_inc)
    elif args.mode == 'allocate':
        write_allocate(fields, args.output_inc)
    elif args.mode == 'deallocate':
        write_deallocate(fields, args.output_inc)
    elif args.mode == 'column_accessor':
        write_column_accessor(fields, args.output_inc)
    elif args.mode == '2d_scalar_accessor':
        write_2d_scalar_accessor(fields, args.output_inc)
    elif args.mode == 'scalar_accessor':
        write_scalar_accessor(fields, args.output_inc)
    elif args.mode == 'column_accessor_int':
        write_column_accessor_int(fields, args.output_inc)
    elif args.mode == '2d_scalar_accessor_int':
        write_2d_scalar_accessor_int(fields, args.output_inc)
    elif args.mode == 'scalar_accessor_int':
        write_scalar_accessor_int(fields, args.output_inc)
    elif args.mode == 'column_accessor_logical':
        write_column_accessor_logical(fields, args.output_inc)
    elif args.mode == '2d_scalar_accessor_logical':
        write_2d_scalar_accessor_logical(fields, args.output_inc)
    elif args.mode == 'scalar_accessor_logical':
        write_scalar_accessor_logical(fields, args.output_inc)
    elif args.mode == 'set_field_2d_real':
        write_set_field_2d_real(fields, args.output_inc)
    elif args.mode == 'set_field_2d_int':
        write_set_field_2d_int(fields, args.output_inc)
    elif args.mode == 'set_field_2d_logical':
        write_set_field_2d_logical(fields, args.output_inc)
    elif args.mode == 'set_field_3d_real':
        write_set_field_3d_real(fields, args.output_inc)
    elif args.mode == 'set_field_3d_int':
        write_set_field_3d_int(fields, args.output_inc)
    elif args.mode == 'set_field_3d_logical':
        write_set_field_3d_logical(fields, args.output_inc)
    elif args.mode == 'set_field_scalar_real':
        write_set_field_scalar_real(fields, args.output_inc)
    elif args.mode == 'set_field_scalar_int':
        write_set_field_scalar_int(fields, args.output_inc)
    elif args.mode == 'set_field_scalar_logical':
        write_set_field_scalar_logical(fields, args.output_inc)
    elif args.mode == 'multiple_fields_interface':
        write_multiple_fields_interface(fields, args.output_inc)
    elif args.mode == 'virtualcolumn_populate':
        write_virtualcolumn_populate(fields, args.output_inc)
    elif args.mode == 'virtualmet_populate':
        write_virtualmet_populate(fields, args.output_inc)
    elif args.mode == 'virtualmet_type':
        write_virtualmet_type(fields, args.output_inc)
    elif args.mode == 'virtualmet_cleanup':
        write_virtualmet_cleanup(fields, args.output_inc)
    else:
        parser.error(f"Unknown mode: {args.mode}")

if __name__ == "__main__":
    main()
