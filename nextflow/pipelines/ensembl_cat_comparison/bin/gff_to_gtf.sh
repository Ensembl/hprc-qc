#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: gff_to_gtf.sh <in.gff3[.gz]> <out.gtf>" >&2
  exit 2
fi
in_gff="$1"
out_gtf="$2"
mkdir -p "$(dirname "$out_gtf")"
# Prepare plain GFF3 for gffread
tmp_in="$(mktemp -t gffread_in.XXXXXX.gff3)"
trap 'rm -f "$tmp_in"' EXIT
case "$in_gff" in
  *.gz) gzip -dc "$in_gff" > "$tmp_in" ;;
  *)    cp -f "$in_gff" "$tmp_in" ;;
 esac
# Run gffread
# -T GTF output; -E tolerate minor issues; -F keep all attrs
gffread -E -F -T -o "$out_gtf" "$tmp_in"
# Lightweight debug log
log="${out_gtf%.gtf}.convert.log"
exons=$(grep -F $'\texon\t' "$out_gtf" | wc -l || true)
tids=$(grep -F 'transcript_id "' "$out_gtf" | wc -l || true)
{
  echo "input: $in_gff"
  echo "gtf: $out_gtf"
  echo "exon_lines=$exons"
  echo "transcript_id_tags=$tids"
} > "$log"
