"""Redraw Figures 1 and 2 of the report from scoring/results (same layout as the figures in the revised draft).
Figure 1: per-case TM-score difference MICA - EModelX(+AF), sorted by resolution. Figure 2: mean TM-score by band."""
import csv, json, statistics as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = "/Users/anony/Downloads/MICA_GAIVI_run"; rounds = {r["case_id"]: r for r in csv.DictReader(open(f"{R}/benchmark/case_rounds.csv"))}
M = ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt"]; res = {m: json.load(open(f"{R}/scoring/results/{m}.json")) for m in M}
tm = lambda m, c: round(res[m][c]["metrics"]["tm_score"], 3)          # values as displayed in Table 1
cases = sorted(rounds, key=lambda c: float(rounds[c]["resolution_A"]))
d = [tm("MICA", c) - tm("EModelX(+AF)", c) for c in cases]; mean = st.mean(d)
fig, ax = plt.subplots(figsize=(8.19, 5.11), dpi=200); y = list(range(len(cases)))[::-1]
for yi, di, c in zip(y, d, cases):
    hi = c.split("_")[1] in ("24VT", "9EFH"); ax.plot([0, di], [yi, yi], color="#9a9a9a", lw=1.2, zorder=1)
    ax.scatter([di], [yi], s=46 if hi else 26, marker="D" if hi else "o", color="#6b6b6b" if hi else "#222b3a", zorder=3)
ax.axvline(0, color="#444444", ls="--", lw=1); ax.axvline(mean, color="#444444", ls=":", lw=1.2)
ax.set_yticks(y); ax.set_yticklabels([f"{c.split('_')[1]}  ({float(rounds[c]['resolution_A']):.2f} Å)" for c in cases], fontsize=7.5)
ax.set_xlabel("TM-score difference: MICA − EModelX(+AF)", fontsize=9); ax.tick_params(axis="x", labelsize=8)
ax.set_xlim(-0.45, 0.45); ax.grid(axis="x", color="#e2e2e2", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{R}/report/figures/image1.png"); plt.close(fig)
bands = ["1.5_to_lt2.5", "2.5_to_lt3.0", "3.0_to_lt3.5", "3.5_to_4.0"]; lab = ["1.5–<2.5", "2.5–<3.0", "3.0–<3.5", "3.5–4.0"]
style = {"MICA": ("-", "o", "#222b3a"), "EModelX(+AF)": ("--", "s", "#555555"), "ModelAngelo": ("-.", "^", "#777777"), "CryoAtom": (":", "D", "#999999"), "EMProt": ((0, (3, 2, 1, 2)), "v", "#4a4a4a")}
fig, ax = plt.subplots(figsize=(7.97, 3.81), dpi=200)
for m in M:
    v = [st.mean(res[m][c]["metrics"]["tm_score"] for c in rounds if rounds[c]["resolution_stratum"] == b) for b in bands]
    ls, mk, col = style[m]; ax.plot(range(4), v, ls=ls, marker=mk, color=col, lw=1.6, ms=6, label=m)
ax.set_xticks(range(4)); ax.set_xticklabels(lab, fontsize=9); ax.set_ylim(0.25, 1.02); ax.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]); ax.tick_params(axis="y", labelsize=9)
ax.set_xlabel("Reported resolution band (Å); six cases per band", fontsize=9); ax.set_ylabel("Mean complex TM-score", fontsize=9)
ax.grid(axis="y", color="#e2e2e2", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(ncol=3, fontsize=8, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(); fig.savefig(f"{R}/report/figures/image2.png"); plt.close(fig)
print("mean difference used for the dotted line: %+.4f" % mean)
