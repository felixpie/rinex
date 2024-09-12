#! /bin/sh
# We want to differentiate 1C and 5Q phase observations
# from two separate receivers, that happen to have
# the same clock.
RINEX_A=OB712480.23O.gz
RINEX_B=gps_10MSps.23O.gz

# We only use the GPS observations to do so.
# We're usually only interested in a single frequency when doing so, but this
# is just an example.
# We reduce to 1 hour dataset
FILTER="GPS;L1C,L5Q"

# Form RINEX(A-B)
# -o: acts as the file prefix
# -q: since we're generating data, we're not interested in opening the workspace
./target/release/rinex-cli \
    -q -o "GPS-L1L5" \
    -P "$FILTER" \
    --fp test_resources/OBS/V3/$RINEX_A \
    diff test_resources/OBS/V3/$RINEX_B

# Export A-B to CSV:
# -o: acts as the CSV file prefix
# -q: since we're generating data, we're not interested in opening the workspace
./target/release/rinex-cli \
    -q -o "GPS-L1L5" \
    -P "$FILTER" \
    --fp test_resources/OBS/V3/$RINEX_A \
    diff --csv test_resources/OBS/V3/$RINEX_B
