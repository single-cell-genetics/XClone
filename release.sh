#!/bin/bash
# XClone Release Script
# Usage: ./release.sh <new_version>
# Example: ./release.sh 0.4.2

set -e  # Exit on error

if [ -z "$1" ]; then
    echo "Usage: ./release.sh <new_version>"
    echo "Example: ./release.sh 0.4.1"
    exit 1
fi

NEW_VERSION=$1
CURRENT_VERSION=$(grep "__version__" xclone/version.py | cut -d'"' -f2)

echo "Current version: $CURRENT_VERSION"
echo "New version: $NEW_VERSION"
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
fi

# Step 1: Update version in version.py
echo "Updating xclone/version.py..."
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" xclone/version.py

# Step 2: Update version in pyproject.toml
echo "Updating pyproject.toml..."
sed -i '' "s/version = \".*\"/version = \"$NEW_VERSION\"/" pyproject.toml

# Step 3: Clean previous builds
echo "Cleaning previous builds..."
rm -rf dist/ build/ *.egg-info

# Step 4: Build distribution
echo "Building distribution..."
python setup.py sdist bdist_wheel

# Step 5: Check version was updated correctly
echo "Verifying version updates..."
if grep -q "__version__ = \"$NEW_VERSION\"" xclone/version.py; then
    echo "✓ version.py updated correctly"
else
    echo "✗ Error: version.py not updated correctly"
    exit 1
fi

if grep -q "version = \"$NEW_VERSION\"" pyproject.toml; then
    echo "✓ pyproject.toml updated correctly"
else
    echo "✗ Error: pyproject.toml not updated correctly"
    exit 1
fi

echo ""
echo "=========================================="
echo "Release $NEW_VERSION prepared successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Update docs/release.rst with release notes"
echo "2. Review changes: git diff"
echo "3. Commit: git commit -am 'Release version $NEW_VERSION'"
echo "4. Tag: git tag -a v$NEW_VERSION -m 'Release version $NEW_VERSION'"
echo "5. Push: git push origin main && git push origin v$NEW_VERSION"
echo "6. Upload to PyPI: twine upload dist/*"
echo ""
