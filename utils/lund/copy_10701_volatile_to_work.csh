#!/bin/tcsh

# copy_10701_volatile_to_work.csh
#
# Copies .hipo simulation output from volatile to work disk.
# Source is NOT deleted.

set SRC = "/lustre24/expphy/volatile/clas12/osg/bulgakov/10701"
set DST = "/w/hallb-scshelf2102/clas12/bulgakov/archive/10701"

echo "Source:"
echo "  $SRC"
echo "Destination:"
echo "  $DST"
echo ""

if ( ! -d "$SRC" ) then
    echo "ERROR: Source directory does not exist:"
    echo "  $SRC"
    exit 1
endif

mkdir -p "$DST"

if ( $status != 0 ) then
    echo "ERROR: Could not create destination directory:"
    echo "  $DST"
    exit 2
endif

echo "Checking source size..."
du -sh "$SRC"

echo ""
echo "Checking destination filesystem space..."
df -h "$DST"

echo ""
echo "Counting source .hipo files..."
set NSRC = `find "$SRC" -type f -name "*.hipo" | wc -l`
echo "Source .hipo files: $NSRC"

echo ""
echo "Starting rsync copy..."
echo "This will COPY files only. It will NOT delete source files."
echo ""

rsync -avh --progress "$SRC"/ "$DST"/

if ( $status != 0 ) then
    echo ""
    echo "ERROR: rsync failed or was interrupted."
    echo "You can rerun this same script; rsync will resume/skip already copied files."
    exit 3
endif

echo ""
echo "Counting destination .hipo files..."
set NDST = `find "$DST" -type f -name "*.hipo" | wc -l`
echo "Destination .hipo files: $NDST"

echo ""
echo "Final sizes:"
echo "Source:"
du -sh "$SRC"
echo "Destination:"
du -sh "$DST"

echo ""
if ( "$NSRC" == "$NDST" ) then
    echo "SUCCESS: File counts match."
else
    echo "WARNING: File counts do not match."
    echo "Source:      $NSRC"
    echo "Destination: $NDST"
endif

echo ""
echo "Original source folder was NOT deleted:"
echo "  $SRC"