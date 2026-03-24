from collections import OrderedDict
from pathlib import Path

infile = "2D_cylinder_beam.data"
outfile = "2D_cylinder_beam_merged.data"

lines = Path(infile).read_text().splitlines()

def find_exact(target):
    for i, l in enumerate(lines):
        if l.strip() == target:
            return i
    raise RuntimeError(f"Cannot find section: {target}")

i_masses = find_exact("Masses")
i_bcoeff = next(i for i, l in enumerate(lines) if l.strip().startswith("Bond Coeffs"))
i_icoeff = next(i for i, l in enumerate(lines) if l.strip().startswith("Improper Coeffs"))
i_atoms  = find_exact("Atoms")
i_bonds  = find_exact("Bonds")
i_impr   = find_exact("Impropers")

# 读取 Bond Coeffs
bond_coeffs = {}
for i in range(i_bcoeff + 2, i_icoeff):
    s = lines[i].strip()
    if not s:
        continue
    p = s.split()
    bond_coeffs[int(p[0])] = (float(p[1]), float(p[2]))   # K, r0

# 读取 Bonds
bonds = []
for i in range(i_bonds + 2, i_impr):
    s = lines[i].strip()
    if not s:
        continue
    p = s.split()
    bonds.append((int(p[0]), int(p[1]), int(p[2]), int(p[3])))  # id, type, a, b

# 合并重复 bond：同一无序原子对 (a,b)
merged = OrderedDict()
for _, btype, a, b in bonds:
    K, r0 = bond_coeffs[btype]
    key = tuple(sorted((a, b)))
    if key not in merged:
        merged[key] = [K, r0, a, b]
    else:
        oldK, oldr0, olda, oldb = merged[key]
        if abs(oldr0 - r0) > 1e-12:
            raise RuntimeError(f"Duplicate bond pair {key} has inconsistent r0: {oldr0} vs {r0}")
        merged[key][0] += K

# 生成新的 Bond Coeffs / Bonds
new_bcoeff = ["Bond Coeffs #harmonic", ""]
new_bonds  = ["Bonds", ""]

for new_id, (_, (K, r0, a, b)) in enumerate(merged.items(), start=1):
    new_bcoeff.append(f"{new_id} {K:.12g} {r0:.12g}")
    new_bonds.append(f"{new_id} {new_id} {a} {b}")

# 更新头部中的 bond 数量
out = lines[:]
for i, l in enumerate(out[:i_masses]):
    s = l.strip()
    if s.endswith(" bonds"):
        out[i] = f"{len(merged)} bonds"
    elif s.endswith(" bond types"):
        out[i] = f"{len(merged)} bond types"

# 重建文件
rebuilt = []
rebuilt.extend(out[:i_bcoeff])
rebuilt.extend(new_bcoeff)
rebuilt.extend(out[i_icoeff:i_bonds])
rebuilt.extend(new_bonds)
rebuilt.extend(out[i_impr:])

Path(outfile).write_text("\n".join(rebuilt) + "\n")

print(f"Written: {outfile}")
print(f"Original bonds: {len(bonds)}")
print(f"Merged bonds:   {len(merged)}")
print(f"Removed duplicates: {len(bonds) - len(merged)}")