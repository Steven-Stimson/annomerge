#!/usr/bin/env python3
"""Poll GitHub device-flow authorization; write token to .gh_token on success."""
import json, time, urllib.request, ssl, sys

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

dev = json.load(open("/work/shidishen/AnnoMerge/annomerge/scripts/.gh_device.json"))
interval = int(dev.get("interval", 5))
url = "https://github.com/login/oauth/access_token"
body = ("grant_type=urn:ietf:params:oauth:grant-type:device_code"
        f"&device_code={dev['device_code']}"
        "&client_id=178c6fc778ccc68e1d6a").encode()

deadline = time.time() + 850
while time.time() < deadline:
    req = urllib.request.Request(url, data=body, headers={
        "Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=30) as r:
            d = json.load(r)
    except Exception as e:
        print("poll error:", e, flush=True)
        time.sleep(interval)
        continue
    if d.get("access_token"):
        with open("/work/shidishen/AnnoMerge/annomerge/scripts/.gh_token", "w") as fh:
            fh.write(d["access_token"])
        print("AUTH_OK", flush=True)
        sys.exit(0)
    err = d.get("error")
    if err == "authorization_pending":
        time.sleep(interval)
    elif err == "slow_down":
        interval += 5
        time.sleep(interval)
    else:
        print("AUTH_FAIL", d, flush=True)
        sys.exit(1)
print("AUTH_TIMEOUT", flush=True)
sys.exit(1)
