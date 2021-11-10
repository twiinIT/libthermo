import json
import matplotlib.pyplot as plt


results = dict()

with open('benchs/benchs_numpy_results.json', 'r') as f:
    results = json.load(f)

fig, (ax1, ax2) = plt.subplots(2, 1)

tail = 10
for y in ("cp", "gamma", "h", "phi", "pr", "eff_poly"):
    ax1.semilogx(results["size"], results[y])
    ax2.semilogx(results["size"][-tail:], results[y][-tail:])

fig.set_size_inches(8, 6)
plt.savefig('benchs/benchs_numpy_results.png', dpi=200)
