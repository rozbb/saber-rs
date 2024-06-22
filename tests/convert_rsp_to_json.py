import json
import re
import sys

usage = """USAGE:
python3 convert_rsp_to_json.py RSPFILE > JSONFILE\
"""

if len(sys.argv) < 2 or "--help" in sys.argv or "-h" in sys.argv:
    print(usage)
    sys.exit(1)

f = open(sys.argv[1], "r")

# Skip the first two lines
next(f)
next(f)

# Accumulate test vectors. These are just dicts of the KAT values
test_vectors = []
while True:
    try:
        test_vectors.append({
            "count": int(re.match(r"count = (\d+)", next(f)).group(1)),
            "seed": re.match(r"seed = (.*)", next(f)).group(1),
            "pk": re.match(r"pk = (.*)", next(f)).group(1),
            "sk": re.match(r"sk = (.*)", next(f)).group(1),
            "ct": re.match(r"ct = (.*)", next(f)).group(1),
            "ss": re.match(r"ss = (.*)", next(f)).group(1),
        })
        # Skip the empty line
        next(f)

    except:
        # Reached end of file
        break

print(json.dumps(test_vectors, indent=4))

f.close()
