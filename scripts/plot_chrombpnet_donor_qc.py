#!/usr/bin/env python3
"""Create a dependency-free SVG summary of donor-resolved ChromBPNet QC."""

import argparse
import json
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def panel(parts, data, x0, title, metric, upper, lower_is_better=False):
    y0, width, height = 55, 450, 280
    parts.extend([
        '<text x="{}" y="26" class="title">{}</text>'.format(x0, title),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#333"/>'.format(x0, y0 + height, x0 + width, y0 + height),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#333"/>'.format(x0, y0, x0, y0 + height),
    ])
    for tick in (0, upper / 2, upper):
        y = y0 + height * (upper - tick) / upper
        parts.extend([
            '<line x1="{}" y1="{:.1f}" x2="{}" y2="{:.1f}" stroke="#ddd"/>'.format(x0, y, x0 + width, y),
            '<text x="{}" y="{:.1f}" text-anchor="end" class="axis">{:.2f}</text>'.format(x0 - 7, y + 4, tick),
        ])
    for i, donor in enumerate(DONORS):
        values = data["donors"][donor]
        model = values["model_metrics"]["peaks"][metric]
        bias = values["matched_scaled_bias_metrics"]["peaks"][metric]
        gx = x0 + 45 + i * 100
        for j, (value, color) in enumerate(((model, "#386cb0"), (bias, "#bdbdbd"))):
            bar_height = height * value / upper
            parts.append('<rect x="{}" y="{:.1f}" width="27" height="{:.1f}" fill="{}"/>'.format(gx + j * 30, y0 + height - bar_height, bar_height, color))
        parts.append('<text x="{}" y="{}" text-anchor="middle" class="axis">{}</text>'.format(gx + 28, y0 + height + 20, donor))
    direction = "lower is better" if lower_is_better else "higher is better"
    parts.append('<text x="{}" y="{}" text-anchor="middle" class="axis">{}</text>'.format(x0 + width / 2, y0 + height + 43, direction))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input_json.read_text())
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1040" height="410" viewBox="0 0 1040 410">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:17px;font-weight:bold}.axis{font-size:12px}.legend{font-size:13px}</style>',
    ]
    panel(parts, data, 65, "Held-out peak count Pearson r", "counts_pearsonr", 0.75)
    panel(parts, data, 570, "Held-out peak profile JSD", "median_jsd", 0.75, lower_is_better=True)
    parts.extend([
        '<rect x="390" y="390" width="14" height="14" fill="#386cb0"/><text x="410" y="402" class="legend">ChromBPNet</text>',
        '<rect x="510" y="390" width="14" height="14" fill="#bdbdbd"/><text x="530" y="402" class="legend">scaled bias</text>',
        '</svg>',
    ])
    args.output_svg.write_text("\n".join(parts) + "\n")


if __name__ == "__main__":
    main()
