#!/bin/bash
#SBATCH --job-name=na2026_slp
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/slp.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --account=mh0033
#
# JJA mean sea-level pressure (hPa) over the North Atlantic-European sector
# (SLP_BOX, 20-80N, 90W-40E) on a 2.5 deg grid, for the NAO / EA EOFs of
# scripts/circulation/01_nao_ea_eof.py and the SLP (blocking) index of
# scripts/heatwave_statistics/fig28_hwd_sst_blocking_regression.py. Absolute JJA
# means, no anomalies: those are formed in the scripts that use them.
#
#   ERA5    monthly MSL (code 151, reduced Gaussian N320), 1940-2026. The
#           archived monthly stream holds 2026 only to July, so 2026 comes from
#           the ERA5T monthly stream entirely (as in 09_era5_zg.sh).
#   MPI-GE  CMIP6 Amon/psl, 50 members, historical + ssp245, 1935-2100: the
#           10-year windows (y-5 to y+4) of 1940 and 2026 need the years around,
#           fig 28 the late century (2061-2100).
#   ICON-EPOC epoc2_020 atm_2d_1mth_mean `psl` on the native R02B08 grid 0054
#           (attached with -setgrid), 1990-2025. Files are stamped at the END of
#           their month, so the months are picked by their run_<YYYYMM> directory.
#   EERIE   r1: CMOR Amon/psl (0.25 deg), hist-1950 + highres-future-ssp245;
#           r2 / r3 (erc2023 / erc2024): the DKRZ 0.25 deg remap of the monthly
#           atm_2d_1mth_mean stream, picked by run_<YYYYMM> directory as well;
#           1975-2030 for all three.
#   MPI-ER  r2 only: ECHAM Amon2d code 151 (Pa) of ER3, hist + ssp585 (the
#           ssp585 files are also named *-hist_*), 1950-2099. Months stamped at
#           their end, so -selmon,6/8 is JJA. ER5 (r3) Amon2d ends 1989 and ER
#           (r1) starts 2014, as for Z300 in 10_model_zg.sh.
#
# Output: data/circulation/<dataset>_jja_slp_2p5deg.nc, variable slp (hPa):
# ERA5 (year, lat, lon), MPI-GE / ICON-EPOC / EERIE / MPI-ER (member, year, lat, lon),
# from stack_slp.py.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

export SLP_BOX=-90,40,20,80
export MPIGE_Y0=1935 MPIGE_Y1=2100
export WORK=${WORK_BASE}/slp
export OUT_SLP=${PROJECT_ROOT}/data/circulation
mkdir -p "${WORK}/era5" "${WORK}/mpige" "$OUT_SLP" "${PROJECT_ROOT}/logs"

echo "=== 1/6 ERA5 JJA MSL, 1940-${REF_YEAR} ==="
era5_year() {
    year=$1
    out=${WORK}/era5/slp_jja_${year}.nc
    [ -s "$out" ] && return 0
    if [ "$year" -eq "$REF_YEAR" ]; then
        src=/pool/data/ERA5/ET/sf/an/1M/151/ETsf00_1M_${year}_151.grb
    else
        src=/pool/data/ERA5/E5/sf/an/1M/151/E5sf00_1M_${year}_151.grb
    fi
    [ "$(cdo -s ntime -selmon,6/8 "$src")" -eq 3 ] || { echo "ERROR: ${src} lacks JJA" >&2; return 1; }
    cdo -s -O -f nc -divc,100 -timmean -sellonlatbox,${SLP_BOX} -remapcon,r144x73 \
        -setgridtype,regular -selmon,6/8 "$src" "$out"
}
export -f era5_year
parallel -j 16 era5_year ::: $(seq 1940 "$REF_YEAR")
cdo -s -O -setattribute,slp@units=hPa -setname,slp -mergetime ${WORK}/era5/slp_jja_*.nc \
    "${WORK}/era5_jja_slp.nc"

echo "=== 2/6 MPI-GE JJA psl, 50 members, ${MPIGE_Y0}-${MPIGE_Y1} ==="
mpige_member() {
    mem=r$1i1p1f1
    out=${WORK}/mpige/slp_jja_${mem}.nc
    [ -s "$out" ] && return 0
    files=""
    for exp in CMIP/MPI-M/MPI-ESM1-2-LR/historical ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245; do
        for f in /pool/data/CMIP6/data/${exp}/${mem}/Amon/psl/gn/v*/psl_*.nc; do
            # keep the files that overlap the years (name ends _YYYYMM-YYYYMM.nc)
            span=${f##*_}; a=${span:0:4}; b=${span:7:4}
            [ "$b" -ge "$MPIGE_Y0" ] && [ "$a" -le "$MPIGE_Y1" ] && files="$files $f"
        done
    done
    # the file holds only Jun/Jul/Aug, so yearmean == JJA mean
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -yearmean \
        -sellonlatbox,${SLP_BOX} -remapcon,r144x73 -selyear,${MPIGE_Y0}/${MPIGE_Y1} \
        -selmon,6/8 -mergetime $files "$out"
}
export -f mpige_member
parallel -j 25 mpige_member ::: $(seq 1 50)

echo "=== 3/6 ICON-EPOC JJA psl, 1990-2025 ==="
export EPOC=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
export GRID_EPOC=/pool/data/ICON/grids/public/mpim/0054/icon_grid_0054_R02B08_G.nc
export WGT_EPOC=${WORK}/wgt_epoc_r144x73.nc
mkdir -p "${WORK}/epoc"
[ -s "$WGT_EPOC" ] || cdo -s gencon,r144x73 -setgrid,"$GRID_EPOC" -selname,psl \
    $(ls ${EPOC}/run_20000601T*/epoc2_020_atm_2d_1mth_mean_*.nc) "$WGT_EPOC"
epoc_year() {
    year=$1
    out=${WORK}/epoc/slp_jja_${year}.nc
    [ -s "$out" ] && return 0
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -sellonlatbox,${SLP_BOX} \
        -remap,r144x73,"$WGT_EPOC" -setgrid,"$GRID_EPOC" -timmean -selname,psl \
        -mergetime ${EPOC}/run_${year}0{6,7,8}01T*/epoc2_020_atm_2d_1mth_mean_*.nc "$out"
}
export -f epoc_year
if [ ! -s "${WORK}/epoc/slp_jja_r1.nc" ]; then
    parallel -j 16 epoc_year ::: $(seq 1990 2025)
    cdo -s -O -mergetime ${WORK}/epoc/slp_jja_[12]*.nc "${WORK}/epoc/slp_jja_r1.nc"
    rm -f ${WORK}/epoc/slp_jja_[12]*.nc
fi

echo "=== 4/6 EERIE r1-r3 JJA psl, 1975-2030 ==="
mkdir -p "${WORK}/eerie"
CMOR=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
if [ ! -s "${WORK}/eerie/slp_jja_r1.nc" ]; then
    r1=$(for y in $(seq 1975 2030); do
        ls ${CMOR}/{hist-1950,highres-future-ssp245}/r1i1p1f1/Amon/psl/gr/v20240618/psl_*_${y}01-${y}12.nc 2>/dev/null || true
    done)
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -yearmean -sellonlatbox,${SLP_BOX} \
        -remapcon,r144x73 -selname,psl -selmon,6/8 -mergetime $r1 "${WORK}/eerie/slp_jja_r1.nc"
fi
for pair in 2:erc2023 3:erc2024; do
    m=${pair%%:*}; exp=${pair#*:}
    [ -s "${WORK}/eerie/slp_jja_r${m}.nc" ] && continue
    B=/work/bm1344/k202193/ICON/${exp}/postprocessing/interpolation
    files=$(for y in $(seq 1975 2030); do for mo in 06 07 08; do
        ls ${B}/{,SSP245/}atm_2d_1mth_mean_remap025/run_${y}${mo}01T*/*.nc 2>/dev/null || true
    done; done)
    # one file per month, stamped at its end: average each JJA by count, then
    # label it with the year of its first month
    for y in $(seq 1975 2030); do
        fy=$(echo "$files" | grep "/run_${y}0[678]01T" || true)
        [ "$(echo "$fy" | grep -c .)" -eq 3 ] || continue
        cdo -s -O -settaxis,${y}-07-16,12:00:00 -timmean -selname,psl -mergetime $fy \
            "${WORK}/eerie/r${m}_${y}.nc"
    done
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -sellonlatbox,${SLP_BOX} \
        -remapcon,r144x73 -mergetime ${WORK}/eerie/r${m}_*.nc "${WORK}/eerie/slp_jja_r${m}.nc"
    rm -f ${WORK}/eerie/r${m}_*.nc
done

echo "=== 5/6 MPI-ER r2 (ER3) JJA SLP, 1950-2099 ==="
mkdir -p "${WORK}/mpier"
if [ ! -s "${WORK}/mpier/slp_jja_r2.nc" ]; then
    ER=/work/uo0122/u241089/MPIESM
    files=$(for y in $(seq 1950 2099); do for m in 06 07 08; do
        ls ${ER}/ER3-{hist,ssp585}/outdata/echam6/Amon2d_ER3-hist_${y}${m}.grb 2>/dev/null || true
    done; done)
    [ "$(echo "$files" | wc -w)" -eq 450 ] || { echo "ERROR: ER3 Amon2d JJA incomplete" >&2; exit 1; }
    # the files hold only Jun/Jul/Aug, so yearmean == JJA mean; -r: ECHAM's
    # absolute time axis (day as %Y%m%d.%f) is not decoded by xarray
    cdo -s -O -r -f nc -setattribute,slp@units=hPa -setname,slp -divc,100 -yearmean \
        -sellonlatbox,${SLP_BOX} -remapcon,r144x73 -selcode,151 -mergetime $files \
        "${WORK}/mpier/slp_jja_r2.nc"
fi

echo "=== 6/6 stack ==="
activate_python
python3 "${SCRIPTS}/pre_process/stack_slp.py"
echo "=== done ==="
