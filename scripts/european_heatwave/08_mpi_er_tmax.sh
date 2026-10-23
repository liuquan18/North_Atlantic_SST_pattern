#!/bin/bash
#SBATCH --job-name=hw_mpier
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_mpier.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --account=mh0033
#
# MPI-ESM1.2-ER (T127 ECHAM6 / TP6M ocean) daily maximum 2 m temperature,
# April-October, Europe 0.25 deg. The daily ECHAM output exists only on tape
# (DKRZ long-term archive, `atlas`/acli; needs a valid ~/.config/acli/session.json):
#
#   r1  ER   ds00007 ER-hist_atm_2014.tar, ds00008 ER-ssp585_atm_<y>.tar 2015-2099
#            -- 1950-2013 was archived as MONTHLY means only (ER-hist_atm_raw.tar),
#               so r1 has no 1991-2020 reference period
#   r2  ER3  ds00009 ER3-hist_atm_<y>.tar 1950-2099 (ssp585 years filed as "hist")
#   r3  ER5  ds00010 ER5-hist_atm_<y>.tar 1950-2099
#
# Each yearly tar is ~60 GB (22 ECHAM streams, one .sz file per month); only the
# daily g3bid stream is used, ~33 MB a month. t2max (code 201) in g3bid is the
# daily maximum: >= the daily mean temp2 everywhere, Europe-land July mean
# Tmax - Tmean = 4.3 K. Steps are stamped 23:56 of the day they cover, so no
# time shift. (The 201 in the monthly `echam` stream is NOT the monthly max of
# these.) The .sz files are szip-compressed GRIB; cdo 2.5 reads them.
#
# Per year: stage (requested for all years up front, so the tape library can
# order the reads), download, extract the Apr-Oct g3bid files, select 201,
# bilinear onto GRID_EU025, then delete the tar. Download ran at ~360 MB/s
# (60 GB in under 3 min); scratch peak is ~60 GB per concurrent year.
#
# Output: ${HW_WORK}/mpier/r<m>/tmax_<year>.nc (degC)
#
# Settings: MPIER_MEMBERS (default "2 3"), MPIER_FIRST_YEAR / MPIER_LAST_YEAR
# (default 1991 / 2026), MPIER_JOBS concurrent years (default 4).

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true
export PATH=/sw/spack-levante/versity-atlas-3.4.5-vrbxoi/bin:$PATH

LTA=/hsm/doku/bm0944
WORK=${HW_WORK}/mpier
TAPE=${WORK}/tape
MEMBERS=${MPIER_MEMBERS:-2 3}
FIRST_YEAR=${MPIER_FIRST_YEAR:-1991}
LAST_YEAR=${MPIER_LAST_YEAR:-2026}
mkdir -p "$TAPE" $(for m in $MEMBERS; do echo "${WORK}/r${m}"; done)

# member, year -> archived tar (empty if the year has no daily output)
tar_path() {
    m=$1; y=$2; base=${LTA}/DKRZ_LTA_944_ds000
    case $m in
        1) if [ "$y" -eq 2014 ]; then echo "${base}07/MPI-M/MPIESM-1-2-ER/ER-hist/ER-hist_atm_2014.tar"
           elif [ "$y" -ge 2015 ]; then echo "${base}08/MPI-M/MPIESM-1-2-ER/ER-ssp585/ER-ssp585_atm_${y}.tar"
           fi ;;
        2) echo "${base}09/MPI-M/MPIESM-1-2-ER/ER3-hist/ER3-hist_atm_${y}.tar" ;;
        3) echo "${base}10/MPI-M/MPIESM-1-2-ER/ER5-hist/ER5-hist_atm_${y}.tar" ;;
    esac
}

is_online() {
    atlas ls --format json "$1" 2>/dev/null | grep -q '"flags": ".o'
}

one_year() {
    m=$1; y=$2; work=$3; tape=$4; grid=$5; months=$6
    out="${work}/r${m}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    src=$(tar_path "$m" "$y")
    [ -n "$src" ] || { echo "  r${m} ${y}: no daily output archived"; return 0; }
    until is_online "$src"; do sleep 60; done

    dir="${tape}/r${m}_${y}"
    rm -rf "$dir"; mkdir -p "$dir"
    # the progress bar goes to stderr: keep it out of the SLURM log unless the download fails
    atlas download "$src" "${dir}/atm.tar" >"${dir}/download.log" 2>&1 \
        || { tail -5 "${dir}/download.log" | tr '\r' '\n' | tail -5; return 1; }
    prefix=$(basename "$src" | sed -E 's/_atm_[0-9]{4}\.tar$//')
    m0=${months%/*}; m1=${months#*/}
    members=$(for mm in $(seq "$m0" "$m1"); do
        printf '%s_echam6_g3bid_%d%02d.sz ' "$prefix" "$y" "$mm"; done)
    tar -xf "${dir}/atm.tar" -C "$dir" $members
    rm -f "${dir}/atm.tar"

    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -remapbil,"$grid" -mergetime \
        $(for f in $members; do echo "-selcode,201 ${dir}/${f}"; done) "${out}.tmp"
    mv "${out}.tmp" "$out"
    rm -rf "$dir"
    echo "  r${m} ${y}: $(cdo -s ntime "$out") days"
}
export -f tar_path is_online one_year
export LTA

# stage every tar still needed in one request
todo=()
for m in $MEMBERS; do
    for y in $(seq "$FIRST_YEAR" "$LAST_YEAR"); do
        [ -s "${WORK}/r${m}/tmax_${y}.nc" ] && continue
        src=$(tar_path "$m" "$y"); [ -n "$src" ] && todo+=("$src")
    done
done
if [ ${#todo[@]} -gt 0 ]; then
    echo "staging ${#todo[@]} tars"
    atlas stage "${todo[@]}"
fi

parallel -j ${MPIER_JOBS:-4} one_year {1} {2} "$WORK" "$TAPE" "$GRID_EU025" "$HW_MONTHS" \
    ::: $MEMBERS ::: $(seq "$FIRST_YEAR" "$LAST_YEAR")
echo "=== done ==="
