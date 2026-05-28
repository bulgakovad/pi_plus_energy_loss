#!/bin/tcsh

set SRC = "/lustre24/expphy/volatile/clas12/bulgakov/lund_files/piplus_electron_lund_multithread_triang_region"
set OUTDIR = "/lustre24/expphy/volatile/clas12/bulgakov/lund_files/archives"
set OUTFILE = "$OUTDIR/piplus_electron_lund_multithread_triang_region.tar.gz"

mkdir -p "$OUTDIR"

tar -czvf "$OUTFILE" \
    -C "$SRC:h" \
    "$SRC:t"

if ( $status != 0 ) then
    echo "ERROR: tar failed"
    exit 1
endif

echo "Archive created:"
ls -lh "$OUTFILE"

echo ""
echo "First files inside archive:"
tar -tzf "$OUTFILE" | head