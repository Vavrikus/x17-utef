#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

QT_DIR="$SCRIPT_DIR/qt6"

echo "=== Setting up Python virtual environment ==="
source /etc/profile.d/30_meta_modules.sh
module add python/3.9.12-gcc-10.2.1-rg2lpmk

python3 -m venv aqt_env
source aqt_env/bin/activate

echo "=== Installing aqtinstall ==="
pip install aqtinstall

echo "=== Determining Latest LTS Version ==="
# Qt6 LTS minor versions follow a predictable 18-month cycle: 6.2, 6.5, 6.8, 6.11, 6.14...
# 1. List all versions
# 2. Convert spaces to newlines
# 3. Filter for LTS minors
# 4. Sort by version and grab the highest patch
ALL_VERSIONS=$(aqt list-qt linux desktop | tr -s ' ' '\n' | grep -E '^6\.')
QT_VERSION=$(echo "$ALL_VERSIONS" | grep -E '^6\.(2|5|8|11|14|17)\.' | sort -V | tail -n 1)

if [ -z "$QT_VERSION" ]; then
    echo "Error: Could not determine the latest Qt LTS version."
    exit 1
fi
echo "Latest LTS version found: $QT_VERSION"

echo "=== Determining Architecture String ==="
# Query the architecture specifically for the found version to handle the gcc_64 vs linux_gcc_64 change
QT_ARCH=$(aqt list-qt linux desktop --arch "$QT_VERSION" | tr -s ' ' '\n' | grep 'gcc_64' | head -n 1)

if [ -z "$QT_ARCH" ]; then
    echo "Error: Could not determine gcc_64 architecture string for Qt $QT_VERSION."
    exit 1
fi
echo "Architecture string found: $QT_ARCH"

echo "=== Downloading Qt $QT_VERSION binaries ==="
python -m aqt install-qt linux desktop "$QT_VERSION" "$QT_ARCH" --outputdir "$QT_DIR"

echo "=== Linking Latest Installation ==="
# Create a permanent 'latest' symlink so other scripts don't need hardcoded version strings
ln -sfn "$QT_DIR/$QT_VERSION/$QT_ARCH" "$QT_DIR/latest"

echo "=== Cleaning up ==="
deactivate
rm -rf aqt_env

echo "Qt6 installation complete at: $QT_DIR/$QT_VERSION/$QT_ARCH"
echo "Symlinked to: $QT_DIR/latest"