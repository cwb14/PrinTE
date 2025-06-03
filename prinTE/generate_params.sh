#!/bin/bash

PARAMS_FILE="params.list"
rm -f "$PARAMS_FILE"

for i_first in {1..10}; do
    fix_first="${i_first}e-8"

    # Convert fix_first to float
    fix_first_num=$(echo "${i_first} * 10^-8" | bc -l)

    # Compute allowed fix_second range
    fix_second_min=$(echo "$fix_first_num + 10^-9" | bc -l)
    fix_second_max=$(echo "$fix_first_num * 2" | bc -l)

    # Available fix_second candidates
    fix_second_candidates=()
    for ((i_second=14; i_second<=140; i_second+=2)); do
        val=$(echo "scale=10; $i_second / 10 * 10^-8" | bc -l)
        fix_second_candidates+=("$val")
    done

    valid_fix_seconds=()
    for val in "${fix_second_candidates[@]}"; do
        if awk -v a="$val" -v min="$fix_second_min" 'BEGIN{if(a > min) exit 0; else exit 1}'; then
            if awk -v a="$val" -v max="$fix_second_max" 'BEGIN{if(a <= max) exit 0; else exit 1}'; then
                valid_fix_seconds+=("$(printf "%.1e\n" "$val" | sed 's/+//;s/0e/1e/')")
            fi
        fi
    done

    if [ ${#valid_fix_seconds[@]} -eq 0 ]; then continue; fi

    for fix_second in "${valid_fix_seconds[@]}"; do
        for k in 0 3 10; do
            for solo_rate in 70 80 95; do
                echo "fix1_${fix_first}_fix2_${fix_second}_k_${k}_solo_${solo_rate}" >> "$PARAMS_FILE"
            done
        done
    done
done

echo "✅ Generated $PARAMS_FILE with $(wc -l < "$PARAMS_FILE") parameter sets"
