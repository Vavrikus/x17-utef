#!/bin/bash

# Extract and source ROOT (ignoring commented lines)
ROOT_SCRIPT=$(grep -v '^\s*#' ~/.bashrc | grep -m 1 "thisroot.sh" | grep -o '\S*thisroot\.sh' | tr -d "\"'")
ROOT_SCRIPT="${ROOT_SCRIPT/#\~/$HOME}"
DEF_FOLDER_JSON='{ 
      "name": "X17", 
      "path": "..", 
      "settings": {
        "cmake.sourceDirectory": "${workspaceFolder}",
        "cmake.cmakePath": "${workspaceFolder}/.vscode/cmake_wrapper.sh"
      }
    }'

if [ -f "$ROOT_SCRIPT" ]; then
    source "$ROOT_SCRIPT"
    echo "Sourcing ROOT successful."
    
    # Check if this is a local build with accessible source code.
    # Assumes the structure where root_src and root_build are in the same parent directory.
    if [ -d "$ROOTSYS/../root_src" ]; then
        ROOT_PARENT_DIR=$(readlink -f "$ROOTSYS/..")
        WORKSPACE_FOLDERS=$DEF_FOLDER_JSON', { "path": "'"$ROOT_PARENT_DIR"'" }'
        
        # Force clangd to index the entire parent directory, resolving the sandbox issue
        cat <<EOF > "$ROOT_PARENT_DIR/.clangd"
CompileFlags:
  CompilationDatabase: root_build
EOF
        echo "Generated .clangd config for ROOT workspace."
    else
        WORKSPACE_FOLDERS=$DEF_FOLDER_JSON
        echo "Note: Local ROOT source directory not found. Jump-to-implementation will be disabled for ROOT." >&2
    fi
else
    WORKSPACE_FOLDERS=$DEF_FOLDER_JSON
    echo "Warning: Could not automatically find thisroot.sh in ~/.bashrc" >&2
fi

# Extract and source Garfield (ignoring commented lines)
GARFIELD_SCRIPT=$(grep -v '^\s*#' ~/.bashrc | grep -m 1 "setupGarfield.sh" | grep -o '\S*setupGarfield\.sh' | tr -d "\"'")
GARFIELD_SCRIPT="${GARFIELD_SCRIPT/#\~/$HOME}"

if [ -f "$GARFIELD_SCRIPT" ]; then
    source "$GARFIELD_SCRIPT"
    echo "Sourcing Garfield successful."
else
    echo "Warning: Could not automatically find setupGarfield.sh in ~/.bashrc" >&2
fi

# Extract and source Geant4 (ignoring commented lines)
GEANT_SCRIPT=$(grep -v '^\s*#' ~/.bashrc | grep -m 1 "geant4.sh" | grep -o '\S*geant4\.sh' | tr -d "\"'")
GEANT_SCRIPT="${GEANT_SCRIPT/#\~/$HOME}"

if [ -f "$GEANT_SCRIPT" ]; then
    source "$GEANT_SCRIPT"
    echo "Sourcing Geant4 successful."
    # Derive Geant4 include path assuming standard <install>/bin/geant4.sh structure
    GEANT_INCLUDE="${GEANT_SCRIPT%/*}/../include/Geant4"
else
    echo "Warning: Could not automatically find geant4.sh in ~/.bashrc" >&2
    GEANT_INCLUDE=""
fi

# Generate VSCode Multi-Root Workspace dynamically
WORKSPACE_FILE="../.vscode/X17.code-workspace"
cat <<EOF > "$WORKSPACE_FILE"
{
  "folders": [
    $WORKSPACE_FOLDERS
  ],
  "settings": {
    "files.associations": {
      "*.shader": "glsl",
      "*.icc": "cpp",
      "plot_drift.C": "cpp",
      "bfield.C": "cpp",
      "cview_update.C": "cpp",
      "geometry.C": "cpp",
      "graph2d_overlay_test.C": "cpp",
      "c_track_xyz_0.C": "cpp",
      "*.inc": "cpp",
      ".clangd": "yaml",
      ".clang-format": "yaml",
      ".clang-tidy": "yaml"
    },
    "C_Cpp.default.includePath": [
      "$(pwd)/../include",
      "$ROOTSYS/include",
      "$GARFIELD_INSTALL/include",
      "$GEANT_INCLUDE"
    ],
    "clangd.arguments": [
        "--background-index",
        "--header-insertion=iwyu",
        "--query-driver=/usr/bin/gcc,/usr/bin/g++"
    ]
  },
  "extensions": {
    "recommendations": [
      "albertopdrf.root-file-viewer",
      "ms-vscode.cpptools",
      "ms-vscode.cmake-tools",
      "llvm-vs-code-extensions.vscode-clangd"
    ]
  }
}
EOF
echo "Workspace file generated: $WORKSPACE_FILE"

# Hand over execution to the real system CMake
command cmake "$@"

echo "Patching compile_commands.json for clangd to use ROOT build directory..."
python3 -c "
import json
import os

try:
    with open('compile_commands.json', 'r') as f:
        data = json.load(f)

    # Automatically derive the install and build paths based on ROOTSYS
    root_install = os.environ.get('ROOTSYS')
    
    if root_install and 'install' in root_install:
        root_build = root_install.replace('install', 'build') # Adjust if your naming differs
        
        if os.path.exists(root_build):
            for entry in data:
                if 'command' in entry:
                    # Force clangd to look at the build headers instead of install headers
                    entry['command'] = entry['command'].replace(f'{root_install}/include', f'{root_build}/include')

            with open('compile_commands.json', 'w') as f:
                json.dump(data, f, indent=2)
            print('Successfully redirected ROOT paths to build directory for clangd.')
        else:
            print('Note: root_build directory not found. Jump-to-implementation will be disabled.')

except Exception as e:
    print(f'Warning: Could not patch compile_commands.json: {e}')
"