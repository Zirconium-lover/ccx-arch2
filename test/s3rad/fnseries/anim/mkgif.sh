#!/usr/bin/env bash
# mkgif.sh <in.mp4> <out.gif> [width=1280] [fps=25] [colours=96]: endless-loop GIF
# with an optimised palette (flat 3b1b-style graphics need few colours; glow
# and photos want more).
in=$1; out=$2; w=${3:-1280}; fps=${4:-25}; nc=${5:-96}
nice ffmpeg -v error -y -i "$in" -vf "fps=$fps,scale=$w:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=$nc:stats_mode=full[p];[b][p]paletteuse=dither=none" -loop 0 "$out"
