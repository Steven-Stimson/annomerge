#!/bin/bash
# 下载 saccharomycetes_odb10：每次失败整体重下（不用 -C -，避免续传拼接损坏）
set -uo pipefail
DEST=/work/shidishen/AnnoMerge/test_data/db/busco
cd "$DEST"
TAR=saccharomycetes_odb10.2024-01-08.tar.gz
URL="https://busco-data.s3.amazonaws.com/lineages/${TAR}"
EXPECT=$(curl --proxy-insecure -k -sI "$URL" | grep -i content-length | tr -dc 0-9)
echo "expected size: $EXPECT"
ok=0
for i in 1 2 3 4 5; do
    rm -f "$TAR"
    curl --proxy-insecure -kL --retry 2 -o "$TAR" "$URL"
    sz=$(stat -c%s "$TAR" 2>/dev/null || echo 0)
    echo "attempt $i size: $sz"
    if [ "$sz" = "$EXPECT" ] && gzip -t "$TAR" 2>/dev/null; then ok=1; echo "download OK"; break; fi
    sleep 3
done
if [ $ok = 1 ]; then
    tar xzf "$TAR" && echo "extracted:" && du -sh saccharomycetes_odb10
else
    echo "DOWNLOAD FAILED after 5 attempts"
fi
