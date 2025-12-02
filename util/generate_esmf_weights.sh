#!/bin/bash
# Generate ESMF Weight Files using ESMF_RegridWeightGen
# This script demonstrates how to pre-generate weight files for optimal CATChem NUOPC performance

set -e

# Configuration
WEIGHTS_DIR="./weights"
INPUT_CONFIG="catchem_input_config.yml"
GRID_DIR="./grids"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
ESMF Weight File Generation Script for CATChem NUOPC

USAGE:
    $0 [COMMAND] [OPTIONS]

COMMANDS:
    setup       Create directories and check ESMF installation
    generate    Generate weight files using ESMF_RegridWeightGen
    template    Create example grid description files
    check       Check generated weight files
    help        Show this help message

OPTIONS:
    -s GRID     Source grid file (SCRIP or ESMF format)
    -d GRID     Destination grid file (SCRIP or ESMF format)
    -m METHOD   Regridding method (bilinear, conservative, nearest, etc.)
    -o OUTPUT   Output weight file name
    -v          Verbose output

EXAMPLES:
    $0 setup                                    # Setup directories
    $0 template                                 # Create example grid files
    $0 generate -s src.nc -d dst.nc -m bilinear -o weights.nc
    $0 check                                    # Validate generated weights

REGRIDDING METHODS:
    bilinear        Bilinear interpolation
    conservative    First-order conservative
    conservative2   Second-order conservative
    neareststod     Nearest neighbor (source to destination)
    nearestdtos     Nearest neighbor (destination to source)

EOF
}

# Check ESMF installation
check_esmf() {
    print_info "Checking ESMF installation..."

    if ! command -v ESMF_RegridWeightGen &> /dev/null; then
        print_error "ESMF_RegridWeightGen not found in PATH"
        print_info "Please ensure ESMF is properly installed and in your PATH"
        print_info "You may need to: module load esmf  or  source esmf_env.sh"
        return 1
    fi

    # Get ESMF version
    local esmf_version
    if command -v ESMF_Info &> /dev/null; then
        esmf_version=$(ESMF_Info | grep "ESMF_VERSION_STRING" | cut -d'"' -f2 2>/dev/null || echo "unknown")
        print_success "Found ESMF version: $esmf_version"
    else
        print_success "ESMF_RegridWeightGen found (version unknown)"
    fi

    return 0
}

# Setup directories
setup_directories() {
    print_info "Setting up directories..."

    for dir in "$WEIGHTS_DIR" "$GRID_DIR"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir"
            print_success "Created directory: $dir"
        else
            print_info "Directory already exists: $dir"
        fi
    done
}

# Create example grid description files
create_templates() {
    print_info "Creating example grid description files..."

    # Create example SCRIP format grid description
    cat > "$GRID_DIR/example_source_grid.py" << 'EOF'
#!/usr/bin/env python3
"""
Example script to create a SCRIP format grid description file
for use with ESMF_RegridWeightGen
"""

import numpy as np
import netCDF4 as nc

def create_latlon_scrip_grid(filename, nlat, nlon, lat_range=(-90, 90), lon_range=(0, 360)):
    """Create a simple lat/lon grid in SCRIP format"""

    # Create coordinate arrays
    lat_centers = np.linspace(lat_range[0], lat_range[1], nlat)
    lon_centers = np.linspace(lon_range[0], lon_range[1], nlon)

    # Create 2D coordinate arrays
    lon_2d, lat_2d = np.meshgrid(lon_centers, lat_centers)

    # Flatten for SCRIP format
    grid_size = nlat * nlon
    grid_center_lat = lat_2d.flatten()
    grid_center_lon = lon_2d.flatten()

    # Calculate grid corners (simple box method)
    dlat = (lat_range[1] - lat_range[0]) / nlat
    dlon = (lon_range[1] - lon_range[0]) / nlon

    # Corner arrays
    grid_corner_lat = np.zeros((grid_size, 4))
    grid_corner_lon = np.zeros((grid_size, 4))

    for i in range(grid_size):
        lat_c = grid_center_lat[i]
        lon_c = grid_center_lon[i]

        # Four corners: SW, SE, NE, NW
        grid_corner_lat[i, :] = [lat_c - dlat/2, lat_c - dlat/2, lat_c + dlat/2, lat_c + dlat/2]
        grid_corner_lon[i, :] = [lon_c - dlon/2, lon_c + dlon/2, lon_c + dlon/2, lon_c - dlon/2]

    # Grid areas (simple calculation)
    grid_area = np.full(grid_size, dlat * dlon * (np.pi/180)**2)

    # Grid mask (all unmasked)
    grid_imask = np.ones(grid_size, dtype=np.int32)

    # Write SCRIP file
    with nc.Dataset(filename, 'w') as f:
        # Dimensions
        f.createDimension('grid_size', grid_size)
        f.createDimension('grid_corners', 4)
        f.createDimension('grid_rank', 2)

        # Variables
        grid_dims = f.createVariable('grid_dims', 'i', ('grid_rank',))
        grid_center_lat_var = f.createVariable('grid_center_lat', 'd', ('grid_size',))
        grid_center_lon_var = f.createVariable('grid_center_lon', 'd', ('grid_size',))
        grid_corner_lat_var = f.createVariable('grid_corner_lat', 'd', ('grid_size', 'grid_corners'))
        grid_corner_lon_var = f.createVariable('grid_corner_lon', 'd', ('grid_size', 'grid_corners'))
        grid_area_var = f.createVariable('grid_area', 'd', ('grid_size',))
        grid_imask_var = f.createVariable('grid_imask', 'i', ('grid_size',))

        # Attributes
        grid_center_lat_var.units = 'degrees'
        grid_center_lon_var.units = 'degrees'
        grid_corner_lat_var.units = 'degrees'
        grid_corner_lon_var.units = 'degrees'
        grid_area_var.units = 'radians^2'

        # Data
        grid_dims[:] = [nlon, nlat]
        grid_center_lat_var[:] = grid_center_lat
        grid_center_lon_var[:] = grid_center_lon
        grid_corner_lat_var[:] = grid_corner_lat
        grid_corner_lon_var[:] = grid_corner_lon
        grid_area_var[:] = grid_area
        grid_imask_var[:] = grid_imask

        # Global attributes
        f.title = f"SCRIP grid: {nlat}x{nlon} lat/lon"
        f.conventions = "SCRIP"

if __name__ == "__main__":
    # Create example grids
    create_latlon_scrip_grid("source_grid_360x180.nc", 180, 360)  # 1-degree global
    create_latlon_scrip_grid("target_grid_720x360.nc", 360, 720)  # 0.5-degree global
    print("Created example SCRIP grid files")
EOF

    chmod +x "$GRID_DIR/example_source_grid.py"
    print_success "Created example grid creation script: $GRID_DIR/example_source_grid.py"

    # Create example weight generation script
    cat > "$GRID_DIR/generate_example_weights.sh" << 'EOF'
#!/bin/bash
# Example weight generation using ESMF_RegridWeightGen

# Generate example grids first
python3 example_source_grid.py

# Generate bilinear weights
ESMF_RegridWeightGen \
    --source source_grid_360x180.nc \
    --destination target_grid_720x360.nc \
    --weight ../weights/example_bilinear_weights.nc \
    --method bilinear \
    --no_log

# Generate conservative weights
ESMF_RegridWeightGen \
    --source source_grid_360x180.nc \
    --destination target_grid_720x360.nc \
    --weight ../weights/example_conservative_weights.nc \
    --method conserve \
    --no_log

echo "Generated example weight files"
EOF

    chmod +x "$GRID_DIR/generate_example_weights.sh"
    print_success "Created example weight generation script: $GRID_DIR/generate_example_weights.sh"
}

# Generate weight files
generate_weights() {
    local source_grid=""
    local dest_grid=""
    local method="bilinear"
    local output=""
    local extra_args=""

    # Parse generation-specific arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -s|--source)
                source_grid="$2"
                shift 2
                ;;
            -d|--destination)
                dest_grid="$2"
                shift 2
                ;;
            -m|--method)
                method="$2"
                shift 2
                ;;
            -o|--output)
                output="$2"
                shift 2
                ;;
            --pole)
                extra_args="$extra_args --pole $2"
                shift 2
                ;;
            --extrap)
                extra_args="$extra_args --extrap_method $2"
                shift 2
                ;;
            --no-log)
                extra_args="$extra_args --no_log"
                shift
                ;;
            *)
                print_warning "Unknown option for generate: $1"
                shift
                ;;
        esac
    done

    # Validate required arguments
    if [ -z "$source_grid" ] || [ -z "$dest_grid" ] || [ -z "$output" ]; then
        print_error "Missing required arguments for weight generation"
        print_info "Usage: $0 generate -s SOURCE_GRID -d DEST_GRID -o OUTPUT_FILE"
        return 1
    fi

    # Check input files
    if [ ! -f "$source_grid" ]; then
        print_error "Source grid file not found: $source_grid"
        return 1
    fi

    if [ ! -f "$dest_grid" ]; then
        print_error "Destination grid file not found: $dest_grid"
        return 1
    fi

    # Create output directory if needed
    local output_dir
    output_dir=$(dirname "$output")
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
        print_success "Created output directory: $output_dir"
    fi

    # Convert method name to ESMF format
    case $method in
        "bilinear")
            esmf_method="bilinear"
            ;;
        "conservative"|"conserve")
            esmf_method="conserve"
            ;;
        "conservative2"|"conserve2")
            esmf_method="conserve"
            extra_args="$extra_args --src_regional --dst_regional"
            ;;
        "nearest"|"neareststod")
            esmf_method="neareststod"
            ;;
        "nearestdtos")
            esmf_method="nearestdtos"
            ;;
        *)
            print_error "Unknown regridding method: $method"
            print_info "Supported methods: bilinear, conservative, conservative2, nearest, nearestdtos"
            return 1
            ;;
    esac

    print_info "Generating weight file..."
    print_info "  Source: $source_grid"
    print_info "  Destination: $dest_grid"
    print_info "  Method: $method ($esmf_method)"
    print_info "  Output: $output"

    # Run ESMF_RegridWeightGen
    if ESMF_RegridWeightGen \
        --source "$source_grid" \
        --destination "$dest_grid" \
        --weight "$output" \
        --method "$esmf_method" \
        $extra_args; then
        print_success "Successfully generated weight file: $output"
    else
        print_error "Failed to generate weight file"
        return 1
    fi
}

# Check generated weight files
check_weights() {
    print_info "Checking generated weight files..."

    if [ ! -d "$WEIGHTS_DIR" ]; then
        print_warning "Weights directory does not exist: $WEIGHTS_DIR"
        return 0
    fi

    local weight_count
    weight_count=$(find "$WEIGHTS_DIR" -name "*.nc" -type f | wc -l)

    if [ "$weight_count" -eq 0 ]; then
        print_warning "No weight files found in $WEIGHTS_DIR"
        return 0
    fi

    print_success "Found $weight_count weight file(s) in $WEIGHTS_DIR"

    # Check each weight file
    for weight_file in "$WEIGHTS_DIR"/*.nc; do
        if [ -f "$weight_file" ]; then
            local filename
            filename=$(basename "$weight_file")
            local filesize
            filesize=$(du -h "$weight_file" | cut -f1)

            print_info "  $filename ($filesize)"

            # Basic validation using ncdump if available
            if command -v ncdump &> /dev/null; then
                if ncdump -h "$weight_file" | grep -q "remap_matrix"; then
                    print_success "    ✓ Valid ESMF weight file format"
                else
                    print_warning "    ? Unknown weight file format"
                fi
            fi
        fi
    done
}

# Parse command line arguments
VERBOSE=false
COMMAND=""

# First pass: extract command
for arg in "$@"; do
    case $arg in
        setup|generate|template|check|help)
            COMMAND="$arg"
            break
            ;;
    esac
done

# Parse remaining arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        setup|generate|template|check|help)
            COMMAND="$1"
            shift
            break
            ;;
        *)
            # For generate command, pass remaining args
            if [ "$COMMAND" = "generate" ]; then
                break
            fi
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Default command
if [ -z "$COMMAND" ]; then
    COMMAND="help"
fi

# Execute command
case $COMMAND in
    setup)
        check_esmf || exit 1
        setup_directories
        ;;
    template)
        setup_directories
        create_templates
        ;;
    generate)
        check_esmf || exit 1
        generate_weights "$@"
        ;;
    check)
        check_weights
        ;;
    help)
        show_help
        ;;
    *)
        print_error "Unknown command: $COMMAND"
        show_help
        exit 1
        ;;
esac
