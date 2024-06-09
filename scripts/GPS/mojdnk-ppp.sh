#!/bin/sh
DATA_DIR=test_resources
CONF=config/survey/ppp_kf.json

./target/release/rinex-cli \
   -f $DATA_DIR/CRNX/V3/MOJN00DNK_R_20201770000_01D_30S_MO.crx.gz \
   -f $DATA_DIR/NAV/V3/MOJN00DNK_R_20201770000_01D_MN.rnx.gz \
   -f $DATA_DIR/SP3/GRG0MGXFIN_20201770000_01D_15M_ORB.SP3.gz \
   -f $DATA_DIR/CLK/V3/GRG0MGXFIN_20201770000_01D_30S_CLK.CLK.gz \
   -P GPS -p -c $CONF | tee logs/mojn-gps+ppp.txt
